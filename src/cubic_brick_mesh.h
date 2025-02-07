//
// Created by Oscar Holroyd on 07/02/2025.
//

#ifndef CUBIC_BRICK_MESH_H
#define CUBIC_BRICK_MESH_H

// Include the OOMPH-LIB header files
#include "generic/brick_mesh.h"
#include "generic/mesh.h"

template <class ELEMENT> class CubicBrickMesh : public virtual BrickMeshBase {
protected:
  /// Number of elements in x direction
  unsigned Nx;

  /// Number of elements in y direction
  unsigned Ny;

  /// Number of elements in z direction
  unsigned Nz;

  /// Number of (linear) points in the element
  unsigned Np;

  /// Minimum value of x coordinate
  double Xmin;

  /// Maximum value of x coordinate
  double Xmax;

  /// Minimum value of y coordinate
  double Ymin;

  /// Minimum value of y coordinate
  double Ymax;

  /// Minimum value of z coordinate
  double Zmin;

  /// Maximum value of z coordinate
  double Zmax;

  /// Boolean to determine if the mesh is periodic in the x direction
  bool Xperiodic;

  /// Boolean to determine if the mesh is periodic in the y direction
  bool Yperiodic;

  /// Generic mesh construction function: contains all the hard work
  void build_mesh(TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper);

  /// Helper function to identify and set the periodically corresponding node
  void set_periodic_node(Node *node_pt, unsigned long el_num,
                         unsigned direction);

public:
  /// Constructor: Pass the number of elements in the x, y, and z directions
  CubicBrickMesh(const unsigned &nx, const unsigned &ny, const unsigned &nz,
                 const double &lx, const double &ly, const double &lz,
                 TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper)
      : Nx(nx), Ny(ny), Nz(nz), Xmin(0.0), Xmax(lx), Ymin(0.0), Ymax(ly),
        Zmin(0.0), Zmax(lz), Xperiodic(false), Yperiodic(false) {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Call the generic build function
    build_mesh(time_stepper_pt);
  }

  /// Constructor: Pass the number of elements in the x, y, and z directions
  CubicBrickMesh(const unsigned &nx, const unsigned &ny, const unsigned &nz,
                 const double &lx, const double &ly, const double &lz,
                 const bool &periodic_in_x, const bool &periodic_in_y,
                 TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper)
      : Nx(nx), Ny(ny), Nz(nz), Xmin(0.0), Xmax(lx), Ymin(0.0), Ymax(ly),
        Zmin(0.0), Zmax(lz), Xperiodic(periodic_in_x),
        Yperiodic(periodic_in_y) {
    // Mesh can only be built with 3D Qelements.
    MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

    // Call the generic build function
    build_mesh(time_stepper_pt);
  }

  /// Access function for number of elements in x directions
  const unsigned &nx() const { return Nx; }

  /// Access function for number of elements in y directions
  const unsigned &ny() const { return Ny; }

  /// Access function for number of elements in y directions
  const unsigned &nz() const { return Nz; }

  /// Return the minimum value of x coordinate
  double x_min() const {
    // Return the value of Xmin
    return Xmin;
  }

  /// Return the maximum value of x coordinate
  double x_max() const {
    // Return the value of Xmax
    return Xmax;
  }

  /// Return the minimum value of y coordinate
  double y_min() const {
    // Return the value of Ymin
    return Ymin;
  }

  /// Return the maximum value of y coordinate
  double y_max() const {
    // Return the value of Ymax
    return Ymax;
  }

  /// Return the minimum value of z coordinate
  double z_min() const {
    // Return the value of Zmin
    return Zmin;
  }

  /// Return the maximum value of z coordinate
  double z_max() const {
    // Return the value of Zmax
    return Zmax;
  }
};

/// given a node on a boundary, the element number containing it and periodic
/// direction (0 for x, 1 for y, 2 for z), this function returns the pointer to
/// the periodic image of the node on the other side of the domain
template <class ELEMENT>
void CubicBrickMesh<ELEMENT>::set_periodic_node(Node *node_pt,
                                                unsigned long el_num,
                                                unsigned direction) {
  // extract the coordinates of the node
  const double x = node_pt->x(0);
  const double y = node_pt->x(1);
  const double z = node_pt->x(2);

  // can only be periodic in the x or y directions
  assert(direction == 0 || direction == 1);

  // the node point must be on the far boundary in the direction specified
  if (direction == 0) {
    // x-periodic, node must be on right boundary (2)
    assert(x == Xmax);
  } else {
    // y-periodic, node must be on rear boundary (3)
    assert(y == Ymax);
  }

  // find the corresponding element and node coordinates on the other side of
  // the domain
  FiniteElement *op_element_pt;
  double op_x, op_y, op_z;
  if (direction == 0) {
    // x-periodic, so go back by one row of elements
    op_element_pt = finite_element_pt(el_num - (Nx - 1));
    op_x = Xmin;
    op_y = y;
    op_z = z;
  } else {
    // y-periodic, so go back by a column of rows of elements
    op_element_pt = finite_element_pt(el_num - (Ny - 1) * Nx);
    op_x = x;
    op_y = Ymin;
    op_z = z;
  }

  // catch the corner case where the element is in the x/y corner
  if (x == Xmax && y == Ymax) {
    op_element_pt = finite_element_pt(el_num - (Ny * Nx - 1));
    op_x = Xmin;
    op_y = Ymin;
    op_z = z;
  }

  // find the node within that element that corresponds to the periodic image
  // of the supplied node
  Node *op_node_pt = nullptr;
  for (unsigned i = 0; i < op_element_pt->nnode(); ++i) {
    Node *n = op_element_pt->node_pt(i);
    if (n->x(0) == op_x && n->x(1) == op_y && n->x(2) == op_z) {
      op_node_pt = n;
      break;
    }
  }

  // ensure that we found the node
  assert(op_node_pt != nullptr);

  // set the periodic node
  node_pt->make_periodic(op_node_pt);
}

//===========================================================================
/// Generic mesh construction. This function contains the "guts" of the
/// mesh generation process, including all the tedious loops, counting
/// and spacing functions. The function should be called in all constuctors
/// of any derived classes.
//===========================================================================
template <class ELEMENT>
void CubicBrickMesh<ELEMENT>::build_mesh(TimeStepper *time_stepper_pt) {
  // Mesh can only be built with 3D Qelements.
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  if ((Nx == 1) || (Ny == 1) || (Nz == 1)) {
    std::ostringstream error_message;
    error_message << "SimpleCubicMesh needs at least two elements in each,\n"
                  << "coordinate direction. You have specified \n"
                  << "Nx=" << Nx << "; Ny=" << Ny << "; Nz=" << Nz << std::endl;
    throw OomphLibError(error_message.str(), OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  // Set the number of boundaries
  set_nboundary(6);

  // Allocate the store for the elements
  Element_pt.resize(Nx * Ny * Nz);
  // Create first element
  unsigned element_num = 0;
  Element_pt[element_num] = new ELEMENT;

  // Read out the number of linear points in the element
  Np = dynamic_cast<ELEMENT *>(finite_element_pt(0))->nnode_1d();

  // Can now allocate the store for the nodes
  Node_pt.resize((1 + (Np - 1) * Nx) * (1 + (Np - 1) * Ny) *
                 (1 + (Np - 1) * Nz));

  // Set up geometrical data
  //------------------------
  unsigned long node_count = 0;

  // Set the length of the elments
  double el_length[3] = {(Xmax - Xmin) / double(Nx), (Ymax - Ymin) / double(Ny),
                         (Zmax - Zmin) / double(Nz)};

  // Storage for local coordinate in element
  Vector<double> s_fraction;

  // Now assign the topology
  // Boundaries are numbered:
  //  0 is at the bottom
  // 1 2 3 4 from the front  proceeding anticlockwise
  // 5 is at the top
  // Pinned value are denoted by an integer value 1
  // Thus if a node is on two boundaries, ORing the values of the
  // boundary conditions will give the most restrictive case (pinning)

  // FIRST ELEMENT (lower left corner at z = 0 plane )
  //----------------------------------

  // Set the corner node
  // Storage for the local node number
  unsigned local_node_num = 0;
  // Allocate memory for the node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer from the element
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmin;
  Node_pt[node_count]->x(1) = Ymin;
  Node_pt[node_count]->x(2) = Zmin;

  // Add the node to boundaries 0, 1, and 4
  add_boundary_node(0, Node_pt[node_count]);
  add_boundary_node(1, Node_pt[node_count]);
  add_boundary_node(4, Node_pt[node_count]);
  // Increment the node number
  node_count++;

  // Loop over the other nodes in the first row in the x direction in
  // the z=0 plane
  for (unsigned l2 = 1; l2 < Np; l2++) {
    // Set the local node number
    local_node_num = l2;

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the local fraction of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundary 0 and 1
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);
    // Increment the node number
    node_count++;
  }

  // Loop over the other node columns in the z = 0 plane
  for (unsigned l1 = 1; l1 < Np; l1++) {
    // Set the local node number
    local_node_num = l1 * Np;

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the first node of the row
    //(with boundaries with 0 and 4)
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundaries 0 and 4
    add_boundary_node(4, Node_pt[node_count]);
    add_boundary_node(0, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Loop over the other nodes in the row
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Set the local node number
      local_node_num = l1 * Np + l2;

      // Allocate the memory for the node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundary 0
      add_boundary_node(0, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }
  }

  //---------------------------------------------------------------------

  // Loop over the other node columns in the z direction
  for (unsigned l3 = 1; l3 < Np; l3++) {
    // Set the local node number
    local_node_num = Np * Np * l3;

    // Allocate memory for the node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

    // Add the node to boundaries 1 and 4
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Loop over the other nodes in the first row in the x direction
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Set the local node number
      local_node_num = l2 + Np * Np * l3;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundary 1
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // Set the local node number
      local_node_num = l1 * Np + Np * Np * l3;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the first node of the row (with boundary 4)
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundary 4
      add_boundary_node(4, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Loop over the other nodes in the row
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Set the local node number
        local_node_num = l2 + l1 * Np + Np * Np * l3;

        // Allocate the memory for the node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // No boundary

        // Increment the node number
        node_count++;
      }
    }
  }

  // END OF FIRST ELEMENT

  //----------------------------------------------------------------------

  // CENTRE OF FIRST ROW OF ELEMENTS (PLANE Z = 0)
  //--------------------------------

  // Now loop over the first row of elements, apart from final element
  for (unsigned j = 1; j < (Nx - 1); j++) {
    // Allocate memory for new element
    element_num = j;
    Element_pt[element_num] = new ELEMENT;

    // We will begin with all the nodes at the plane z = 0

    // Do first row of nodes

    // First column of nodes is same as neighbouring element
    finite_element_pt(element_num)->node_pt(0) =
        finite_element_pt(element_num - 1)->node_pt((Np - 1));

    // New nodes for other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Set the local node number
      local_node_num = l2;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 0 an 1
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Do the rest of the nodes at the plane z = 0
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(element_num)->node_pt(l1 * Np) =
          finite_element_pt(element_num - 1)->node_pt(l1 * Np + (Np - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Set the local node number
        local_node_num = l2 + l1 * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the Boundary
        add_boundary_node(0, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }
    }

    // Loop over the other node columns in the z direction
    for (unsigned l3 = 1; l3 < Np; l3++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(j)->node_pt(l3 * Np * Np) =
          finite_element_pt(j - 1)->node_pt(l3 * Np * Np + (Np - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Set the local node number
        local_node_num = l2 + l3 * Np * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j)->node_pt(l1 * Np + l3 * Np * Np) =
            finite_element_pt(j - 1)->node_pt(l1 * Np + (Np - 1) +
                                              l3 * Np * Np);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Set the local node number
          local_node_num = l2 + l1 * Np + l3 * Np * Np;

          // Allocate memory for the nodes
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }
    }
  }

  // MY FINAL ELEMENT IN FIRST ROW (lower right corner)
  //-----------------------------------------------

  // Allocate memory for new element
  element_num = Nx - 1;
  Element_pt[element_num] = new ELEMENT;

  // We will begin with all the nodes at the plane z = 0

  // Do first row of nodes

  // First node is same as neighbouring element
  finite_element_pt(element_num)->node_pt(0) =
      finite_element_pt(element_num - 1)->node_pt((Np - 1));

  // New nodes for other columns
  for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
    // Set the local node number
    local_node_num = l2;

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundaries 0 an 1
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(1, Node_pt[node_count]);
    // Increment the node number
    node_count++;
  }

  // Last node (corner)
  // Set the local node number
  local_node_num = Np - 1;

  // Allocate memory for the node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer from the element
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmax;
  Node_pt[node_count]->x(1) = Ymin;
  Node_pt[node_count]->x(2) = Zmin;

  // if required, make it periodic with the node on the other side
  if (Xperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 0);
  }

  // Add the node to the boundaries
  add_boundary_node(0, Node_pt[node_count]);
  add_boundary_node(1, Node_pt[node_count]);
  add_boundary_node(2, Node_pt[node_count]);
  // Increment the node number
  node_count++;

  // Do the middle nodes at the plane z = 0
  for (unsigned l1 = 1; l1 < Np; l1++) {
    // First column of nodes is same as neighbouring element
    finite_element_pt(element_num)->node_pt(l1 * Np) =
        finite_element_pt(element_num - 1)->node_pt(l1 * Np + (Np - 1));

    // New nodes for other columns
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Set the local node number
      local_node_num = l2 + l1 * Np;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundary 0
      add_boundary_node(0, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // New node for final column
    // Set the local node number
    local_node_num = l1 * Np + (Np - 1);

    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
    Node_pt[node_count]->x(2) = Zmin;

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    }

    // Add the node to the boundaries 0 and 2
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(0, Node_pt[node_count]);
    // Increment the node number
    node_count++;
  }

  // Loop over the other node columns in the z direction
  for (unsigned l3 = 1; l3 < Np; l3++) {
    // First column of nodes is same as neighbouring element
    finite_element_pt(element_num)->node_pt(l3 * Np * Np) =
        finite_element_pt(element_num - 1)->node_pt(l3 * Np * Np + (Np - 1));

    // New nodes for other rows (y=0)
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Set the local node number
      local_node_num = l2 + l3 * Np * Np;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundary 1
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Last node of the row (y=0)
    // Set the local node number
    local_node_num = Np - 1 + l3 * Np * Np;

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    }

    // Add the node to the boundary 1 and 2
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Do the rest of the nodes
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(element_num)->node_pt(l1 * Np + l3 * Np * Np) =
          finite_element_pt(element_num - 1)
              ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Set the local node number
        local_node_num = l2 + l1 * Np + l3 * Np * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // No boundaries

        // Increment the node number
        node_count++;
      }

      // Last nodes (at the surface x = Lx)
      // Set the local nod number
      local_node_num = l1 * Np + (Np - 1) + l3 * Np * Np;
      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundary 2
      add_boundary_node(2, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }
  }

  // ALL CENTRAL ELEMENT ROWS (WE ARE STILL IN THE LAYER z=0)
  //------------------------

  // Loop over remaining element rows
  for (unsigned i = 1; i < (Ny - 1); i++) {
    // Set the first element in the row

    // Allocate memory for element (z = 0) (x =0)
    element_num = Nx * i;
    Element_pt[element_num] = new ELEMENT;

    // The first row of nodes is copied from the element below (at z=0)
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(element_num)->node_pt(l2) =
          finite_element_pt(element_num - Nx)->node_pt((Np - 1) * Np + l2);
    }

    // Other rows are new nodes
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column of nodes
      // Set the local node number
      local_node_num = l1 * Np;

      // Allocate memory for the fist node in the x direction (x=0)
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundaries 4 and 0
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Set the local node number
        local_node_num = l2 + l1 * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary  and 0
        add_boundary_node(0, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }
    }

    // As always we extend now this element to the third direction
    for (unsigned l3 = 1; l3 < Np; l3++) {
      // The first row of nodes is copied from the element below (at z=0)
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(element_num)->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(element_num - Nx)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      // Other rows are new nodes (first the nodes for which x=0)
      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column of nodes
        // Set the local node number
        local_node_num = l1 * Np + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns (we extend to the rest of the nodes)
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Set the local node number
          local_node_num = l2 + l1 * Np + Np * Np * l3;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }
    }

    // Now loop over the rest of the elements in the row, apart from the last
    for (unsigned j = 1; j < (Nx - 1); j++) {
      // Allocate memory for new element
      element_num = Nx * i + j;
      Element_pt[element_num] = new ELEMENT;

      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(element_num)->node_pt(l2) =
            finite_element_pt(element_num - Nx)->node_pt((Np - 1) * Np + l2);
      }

      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column is same as neighbouring element
        finite_element_pt(element_num)->node_pt(l1 * Np) =
            finite_element_pt(element_num - 1)->node_pt(l1 * Np + (Np - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Set local node number
          local_node_num = l1 * Np + l2;

          // Allocate memory for the nodes
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin;

          // Add the node to the boundary  and 0
          add_boundary_node(0, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }
      }

      // We extend to the third dimension
      for (unsigned l3 = 1; l3 < Np; l3++) {
        // The first row is copied from the lower element

        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(element_num)->node_pt(l2 + l3 * Np * Np) =
              finite_element_pt(element_num - Nx)
                  ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
        }

        for (unsigned l1 = 1; l1 < Np; l1++) {
          // First column is same as neighbouring element
          finite_element_pt(element_num)->node_pt(l1 * Np + l3 * Np * Np) =
              finite_element_pt(element_num - 1)
                  ->node_pt(l1 * Np + l3 * Np * Np + (Np - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < Np; l2++) {
            // Set the local node number
            local_node_num = l1 * Np + l2 + l3 * Np * Np;

            // Allocate memory for the nodes
            Node_pt[node_count] =
                finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }

    } // End of loop over elements in row

    // Do final element in row

    // Allocate memory for element
    element_num = Nx * i + Nx - 1;
    Element_pt[element_num] = new ELEMENT;

    // We begin with z =0
    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(element_num)->node_pt(l2) =
          finite_element_pt(element_num - Nx)->node_pt((Np - 1) * Np + l2);
    }

    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column is same as neighbouring element
      finite_element_pt(element_num)->node_pt(l1 * Np) =
          finite_element_pt(element_num - 1)->node_pt(l1 * Np + (Np - 1));

      // Middle nodes
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Set local node number
        local_node_num = l1 * Np + l2;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary  and 0
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Final node

      // Set the local node number
      local_node_num = l1 * Np + (Np - 1);

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin;

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundaries - and 1
      add_boundary_node(0, Node_pt[node_count]);
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;

    } // End of loop over rows of nodes in the element

    // We go to the third dimension
    for (unsigned l3 = 1; l3 < Np; l3++) {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(element_num)->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(element_num - Nx)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column is same as neighbouring element
        finite_element_pt(element_num)->node_pt(l1 * Np + l3 * Np * Np) =
            finite_element_pt(element_num - 1)
                ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
          // Set the local node number
          local_node_num = l1 * Np + l2 + l3 * Np * Np;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Final node
        // Set the local node number
        local_node_num = l1 * Np + (Np - 1) + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // if required, make it periodic with the node on the other side
        if (Xperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 0);
        }

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      } // End of loop over rows of nodes in the element

    } // End of the 3dimension loop

  } // End of loop over rows of elements

  // FINAL ELEMENT ROW (IN THE z=0 LAYER)
  //=================

  // FIRST ELEMENT IN UPPER ROW (upper left corner)
  //----------------------------------------------

  // Allocate memory for element
  element_num = Nx * (Ny - 1);
  Element_pt[element_num] = new ELEMENT;
  // We begin with all the nodes for which z=0
  // The first row of nodes is copied from the element below
  for (unsigned l2 = 0; l2 < Np; l2++) {
    finite_element_pt(element_num)->node_pt(l2) =
        finite_element_pt(element_num - Nx)->node_pt((Np - 1) * Np + l2);
  }

  // Second row of  nodes
  // First column of nodes
  for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
    // Set the local node number
    local_node_num = Np * l1;

    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
    Node_pt[node_count]->x(2) = Zmin;

    // Add the node to the boundaries 4 and 0
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Set the local node number
      local_node_num = Np * l1 + l2;

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundary 0
      add_boundary_node(0, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }
  }

  // Final row of nodes
  // First column of nodes
  // Top left node
  // Set local node number
  local_node_num = Np * (Np - 1);
  // Allocate memory for node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer from the element
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmin;
  Node_pt[node_count]->x(1) = Ymax;
  Node_pt[node_count]->x(2) = Zmin;

  if (Yperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 1);
  }

  // Add the node to the boundaries 0,3 and 4
  add_boundary_node(0, Node_pt[node_count]);
  add_boundary_node(3, Node_pt[node_count]);
  add_boundary_node(4, Node_pt[node_count]);

  // Increment the node number
  node_count++;

  // Now do the other columns
  for (unsigned l2 = 1; l2 < Np; l2++) {
    // Set the local node number
    local_node_num = Np * (Np - 1) + l2;
    // Allocate memory for the node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin;

    if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);
    // Increment the node number
    node_count++;
  }

  // We jump to the third dimension
  for (unsigned l3 = 1; l3 < Np; l3++) {
    // The first row of nodes is copied from the element below
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(element_num)->node_pt(l2 + l3 * Np * Np) =
          finite_element_pt(element_num - Nx)
              ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
    }

    // Second row of  nodes
    // First column of nodes
    for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
      // Set the local node number
      local_node_num = Np * l1 + l3 * Np * Np;

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // Add the node to the boundary 4
      add_boundary_node(4, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Set local node number
        local_node_num = Np * l1 + l2 + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // No boundaries

        // Increment the node number
        node_count++;
      }
    }

    // Final row of nodes
    // First column of nodes
    // Top left node
    local_node_num = Np * (Np - 1) + l3 * Np * Np;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

    if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
      // Allocate memory for the node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundary 3
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }
  }

  // Now loop over the rest of the elements in the row, apart from the last
  // (first plane z = 0)
  for (unsigned j = 1; j < (Nx - 1); j++) {
    // Allocate memory for element
    element_num = Nx * (Ny - 1) + j;
    Element_pt[element_num] = new ELEMENT;
    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(element_num)->node_pt(l2) =
          finite_element_pt(element_num - Nx)->node_pt((Np - 1) * Np + l2);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
      // First column is same as neighbouring element
      finite_element_pt(element_num)->node_pt(Np * l1) =
          finite_element_pt(element_num - 1)->node_pt(Np * l1 + (Np - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        local_node_num = Np * l1 + l2;
        // Allocate memory for the node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);

        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin;

        // Add the node to the boundary
        add_boundary_node(0, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // Top row
    // First column is same as neighbouring element
    finite_element_pt(element_num)->node_pt(Np * (Np - 1)) =
        finite_element_pt(element_num - 1)->node_pt(Np * (Np - 1) + (Np - 1));
    // New nodes for other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      local_node_num = Np * (Np - 1) + l2;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin;

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundary
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Jump in the third dimension

    for (unsigned l3 = 1; l3 < Np; l3++) {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(element_num)->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(element_num - Nx)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
        // First column is same as neighbouring element
        finite_element_pt(element_num)->node_pt(Np * l1 + l3 * Np * Np) =
            finite_element_pt(element_num - 1)
                ->node_pt(Np * l1 + (Np - 1) + l3 * Np * Np);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          local_node_num = Np * l1 + l2 + l3 * Np * Np;
          // Allocate memory for the node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);

          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

          // No boundaries

          // Increment the node number
          // add_boundary_node(0,Node_pt[node_count]);
          node_count++;
        }
      }

      // Top row
      // First column is same as neighbouring element
      finite_element_pt(element_num)->node_pt(Np * (Np - 1) + l3 * Np * Np) =
          finite_element_pt(element_num - 1)
              ->node_pt(Np * (Np - 1) + (Np - 1) + l3 * Np * Np);
      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        if (Yperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 1);
        }

        // Add the node to the boundary
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

  } // End of loop over central elements in row

  // FINAL ELEMENT IN ROW (upper right corner) IN LAYER z = 0
  //--------------------------------------------------------

  // Allocate memory for element
  element_num = Nx * (Ny - 1) + Nx - 1;
  Element_pt[element_num] = new ELEMENT;

  // We work first in the plane z =0
  // The first row is copied from the lower element
  for (unsigned l2 = 0; l2 < Np; l2++) {
    finite_element_pt(element_num)->node_pt(l2) =
        finite_element_pt(element_num - Nx)->node_pt((Np - 1) * Np + l2);
  }

  // Second rows
  for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
    // First column is same as neighbouring element
    finite_element_pt(element_num)->node_pt(Np * l1) =
        finite_element_pt(element_num - 1)->node_pt(Np * l1 + (Np - 1));

    // Middle nodes
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      local_node_num = Np * l1 + l2;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin;

      // Add the node to the boundary
      add_boundary_node(0, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Final node
    local_node_num = Np * l1 + (Np - 1);
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
    Node_pt[node_count]->x(2) = Zmin;

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    }

    // Add the node to the boundaries 0 and 2
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);

    // Increment the node number
    node_count++;

  } // End of loop over middle rows

  // Final row
  // First column is same as neighbouring element
  finite_element_pt(element_num)->node_pt(Np * (Np - 1)) =
      finite_element_pt(element_num - 1)->node_pt(Np * (Np - 1) + (Np - 1));

  // Middle nodes
  for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
    local_node_num = Np * (Np - 1) + l2;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin;

    if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries 0 nd 3
    add_boundary_node(0, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // Final node
  // Determine number of values
  local_node_num = Np * (Np - 1) + (Np - 1);
  // Allocate memory for node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmax;
  Node_pt[node_count]->x(1) = Ymax;
  Node_pt[node_count]->x(2) = Zmin;

  // if required, make it periodic with the node on the other side
  if (Xperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 0);
  } else if (Yperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 1);
  }

  // Add the node to the boundaries 0,2 and 3
  add_boundary_node(0, Node_pt[node_count]);
  add_boundary_node(2, Node_pt[node_count]);
  add_boundary_node(3, Node_pt[node_count]);

  // Increment the node number
  node_count++;

  // We jump to the third dimension

  for (unsigned l3 = 1; l3 < Np; l3++) {
    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(element_num)->node_pt(l2 + l3 * Np * Np) =
          finite_element_pt(element_num - Nx)
              ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
      // First column is same as neighbouring element
      finite_element_pt(element_num)->node_pt(Np * l1 + l3 * Np * Np) =
          finite_element_pt(element_num - 1)
              ->node_pt(Np * l1 + (Np - 1) + l3 * Np * Np);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Determine number of values
        local_node_num = Np * l1 + l2 + l3 * Np * Np;
        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

        // No boundaries

        // Increment the node number
        node_count++;
      }

      // Final node
      // Determine number of values
      local_node_num = Np * l1 + (Np - 1) + l3 * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundary 2
      add_boundary_node(2, Node_pt[node_count]);
      // Increment the node number
      node_count++;

    } // End of loop over middle rows

    // Final row
    // First column is same as neighbouring element
    finite_element_pt(element_num)->node_pt(Np * (Np - 1) + l3 * Np * Np) =
        finite_element_pt(element_num - 1)
            ->node_pt(Np * (Np - 1) + (Np - 1) + l3 * Np * Np);

    // Middle nodes
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Determine number of values
      local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundary 3
      add_boundary_node(3, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Final node
    // Determine number of values
    local_node_num = Np * (Np - 1) + (Np - 1) + l3 * Np * Np;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * s_fraction[2];

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    } else if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries 2 and 3
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // END OF THE FIRST LAYER

  //----------------------------------------------------------------------------------------------------------------------------
  // ***************************************NOW WE MAKE ALL THE INTREMEDIATE
  // LAYERS**********************************************
  //----------------------------------------------------------------------------------------------------------------------------

  for (unsigned k = 1; k < (Nz - 1); k++) // good loop for the diferent layers
  // for(unsigned k=1;k<Nz;k++)  // bad loop for testing the hole mesh but the
  // last layer
  {
    // FIRST ELEMENT OF THE LAYER
    //----------------------------------

    element_num = k * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(element_num)->node_pt(l2 + Np * l1) =
            finite_element_pt(element_num - Nx * Ny)
                ->node_pt(l2 + Np * l1 + Np * Np * (Np - 1));
      }
    }

    // Loop over the other node columns in the z direction

    for (unsigned l3 = 1; l3 < Np; l3++) {
      // Set the corner node
      // Determine number of values at this node
      local_node_num = Np * Np * l3;

      // Allocate memory for the node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);

      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

      // Add the node to boundaries 1 and 4
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Loop over the other nodes in the first row in the x direction
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine the number of values at this node
        local_node_num = l2 + Np * Np * l3;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Loop over the other node columns
      for (unsigned l1 = 1; l1 < Np; l1++) {
        // Determine the number of values
        local_node_num = l1 * Np + Np * Np * l3;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the first node of the row (with boundary 4)
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);
        // Increment the node number
        node_count++;

        // Loop over the other nodes in the row
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Set the number of values
          local_node_num = l1 * Np + l2 + Np * Np * l3;

          // Allocate the memory for the node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // No boundary

          // Increment the node number
          node_count++;
        }
      }
    }

    //----------------------------------------------------------------------

    // CENTRE OF FIRST ROW OF ELEMENTS
    //--------------------------------

    // Now loop over the first row of elements, apart from final element
    for (unsigned j = 1; j < (Nx - 1); j++) {
      // Allocate memory for new element
      element_num = j + k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < Np; l1++) {
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(j + k * Nx * Ny)->node_pt(l2 + Np * l1) =
              finite_element_pt(j + (k - 1) * Nx * Ny)
                  ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
        }
      }

      // Loop over the other node columns in the z direction
      for (unsigned l3 = 1; l3 < Np; l3++) {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j + k * Nx * Ny)->node_pt(l3 * Np * Np) =
            finite_element_pt(j - 1 + k * Nx * Ny)
                ->node_pt(l3 * Np * Np + (Np - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = l2 + l3 * Np * Np;

          // Allocate memory for the nodes
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 1
          add_boundary_node(1, Node_pt[node_count]);
          // Increment the node number
          node_count++;
        }

        // Do the rest of the nodes
        for (unsigned l1 = 1; l1 < Np; l1++) {
          // First column of nodes is same as neighbouring element
          finite_element_pt(j + k * Nx * Ny)->node_pt(l1 * Np + l3 * Np * Np) =
              finite_element_pt(j - 1 + k * Nx * Ny)
                  ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < Np; l2++) {
            // Determine number of values
            local_node_num = l1 * Np + l2 + l3 * Np * Np;

            // Allocate memory for the nodes
            Node_pt[node_count] =
                finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
            Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }
    }

    // MY FINAL ELEMENT IN FIRST ROW (right corner)
    //-----------------------------------------------

    // Allocate memory for new element
    element_num = Nx - 1 + k * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx - 1 + k * Nx * Ny)->node_pt(l2 + Np * l1) =
            finite_element_pt(Nx - 1 + (k - 1) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // Loop over the other node columns in the z direction
    for (unsigned l3 = 1; l3 < Np; l3++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(Nx - 1 + k * Nx * Ny)->node_pt(l3 * Np * Np) =
          finite_element_pt(Nx - 2 + k * Nx * Ny)
              ->node_pt(l3 * Np * Np + (Np - 1));

      // New nodes for other rows (y=0)
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Determine number of values
        local_node_num = l2 + l3 * Np * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Last node of the row (y=0)

      // Determine number of values
      local_node_num = (Np - 1) + l3 * Np * Np;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundary 1 and 2
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(2, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column of nodes is same as neighbouring element
        finite_element_pt(Nx - 1 + k * Nx * Ny)
            ->node_pt(l1 * Np + l3 * Np * Np) =
            finite_element_pt(Nx - 2 + k * Nx * Ny)
                ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
          // Determine number of values
          local_node_num = l1 * Np + l2 + l3 * Np * Np;

          // Allocate memory for the nodes
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Last nodes (at the surface x = Lx)
        // Determine number of values
        local_node_num = l1 * Np + (Np - 1) + l3 * Np * Np;
        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // if required, make it periodic with the node on the other side
        if (Xperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 0);
        }

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // ALL CENTRAL ELEMENT ROWS
    //------------------------

    // Loop over remaining element rows
    for (unsigned i = 1; i < (Ny - 1); i++) {
      // Set the first element in the row

      // Allocate memory for element (x =0)
      element_num = Nx * i + Nx * Ny * k;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < Np; l1++) {
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * i + k * Nx * Ny)->node_pt(l2 + Np * l1) =
              finite_element_pt(Nx * i + (k - 1) * Nx * Ny)
                  ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
        }
      }

      // We extend now this element to the third direction

      for (unsigned l3 = 1; l3 < Np; l3++) {
        // The first row of nodes is copied from the element below (at z=0)
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * i + k * Nx * Ny)->node_pt(l2 + l3 * Np * Np) =
              finite_element_pt(Nx * (i - 1) + k * Nx * Ny)
                  ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
        }

        // Other rows are new nodes (first the nodes for which x=0)
        for (unsigned l1 = 1; l1 < Np; l1++) {
          // First column of nodes

          // Determine number of values
          local_node_num = l1 * Np + l3 * Np * Np;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // Add the node to the boundary 4
          add_boundary_node(4, Node_pt[node_count]);

          // Increment the node number
          node_count++;

          // Now do the other columns (we extend to the rest of the nodes)
          for (unsigned l2 = 1; l2 < Np; l2++) {
            // Determine number of values
            local_node_num = l1 * Np + l2 + Np * Np * l3;

            // Allocate memory for node
            Node_pt[node_count] =
                finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer from the element
            finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
            Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }

      // Now loop over the rest of the elements in the row, apart from the
      // last
      for (unsigned j = 1; j < (Nx - 1); j++) {
        // Allocate memory for new element
        element_num = Nx * i + j + k * Nx * Ny;
        Element_pt[element_num] = new ELEMENT;

        // The lowest nodes are copied from the lower element
        for (unsigned l1 = 0; l1 < Np; l1++) {
          for (unsigned l2 = 0; l2 < Np; l2++) {
            finite_element_pt(Nx * i + j + k * Nx * Ny)->node_pt(l2 + Np * l1) =
                finite_element_pt(Nx * i + j + (k - 1) * Nx * Ny)
                    ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
          }
        }

        // We extend to the third dimension

        for (unsigned l3 = 1; l3 < Np; l3++) {
          // The first row is copied from the lower element

          for (unsigned l2 = 0; l2 < Np; l2++) {
            finite_element_pt(Nx * i + j + k * Nx * Ny)
                ->node_pt(l2 + l3 * Np * Np) =
                finite_element_pt(Nx * (i - 1) + j + k * Nx * Ny)
                    ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
          }

          for (unsigned l1 = 1; l1 < Np; l1++) {
            // First column is same as neighbouring element
            finite_element_pt(Nx * i + j + k * Nx * Ny)
                ->node_pt(l1 * Np + l3 * Np * Np) =
                finite_element_pt(Nx * i + (j - 1) + k * Nx * Ny)
                    ->node_pt(l1 * Np + l3 * Np * Np + (Np - 1));

            // New nodes for other columns
            for (unsigned l2 = 1; l2 < Np; l2++) {
              // Determine number of values
              local_node_num = l1 * Np + l2 + l3 * Np * Np;

              // Allocate memory for the nodes
              Node_pt[node_count] =
                  finite_element_pt(element_num)
                      ->construct_node(local_node_num, time_stepper_pt);
              // Set the pointer
              finite_element_pt(element_num)->node_pt(local_node_num) =
                  Node_pt[node_count];
              // Get the fractional position of the node
              finite_element_pt(element_num)
                  ->local_fraction_of_node(local_node_num, s_fraction);

              // Set the position of the node
              Node_pt[node_count]->x(0) =
                  Xmin + el_length[0] * (j + s_fraction[0]);
              Node_pt[node_count]->x(1) =
                  Ymin + el_length[1] * (i + s_fraction[1]);
              Node_pt[node_count]->x(2) =
                  Zmin + el_length[2] * (k + s_fraction[2]);

              // No boundaries
              // Increment the node number
              node_count++;
            }
          }
        }

      } // End of loop over elements in row

      // Do final element in row

      // Allocate memory for element
      element_num = Nx * i + Nx - 1 + k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < Np; l1++) {
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * i + Nx - 1 + k * Nx * Ny)
              ->node_pt(l2 + Np * l1) =
              finite_element_pt(Nx * i + Nx - 1 + (k - 1) * Nx * Ny)
                  ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
        }
      }

      // We go to the third dimension
      for (unsigned l3 = 1; l3 < Np; l3++) {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * i + Nx - 1 + k * Nx * Ny)
              ->node_pt(l2 + l3 * Np * Np) =
              finite_element_pt(Nx * (i - 1) + Nx - 1 + k * Nx * Ny)
                  ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
        }

        for (unsigned l1 = 1; l1 < Np; l1++) {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + Nx - 1 + k * Nx * Ny)
              ->node_pt(l1 * Np + l3 * Np * Np) =
              finite_element_pt(Nx * i + Nx - 2 + k * Nx * Ny)
                  ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

          // Middle nodes
          for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
            // Determine number of values
            local_node_num = l1 * Np + l2 + l3 * Np * Np;

            // Allocate memory for node
            Node_pt[node_count] =
                finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
            Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            node_count++;
          }

          // Final node

          // Determine number of values
          local_node_num = l1 * Np + (Np - 1) + l3 * Np * Np;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmax;
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // if required, make it periodic with the node on the other side
          if (Xperiodic) {
            set_periodic_node(Node_pt[node_count], element_num, 0);
          }

          // Add the node to the boundary 2
          add_boundary_node(2, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        } // End of loop over rows of nodes in the element

      } // End of the 3dimension loop

    } // End of loop over rows of elements

    // FINAL ELEMENT ROW (IN INTERMIDIATE  LAYERS)
    //=================

    // FIRST ELEMENT IN UPPER ROW (upper left corner)
    //----------------------------------------------

    // Allocate memory for element
    element_num = Nx * (Ny - 1) + k * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * (Ny - 1) + k * Nx * Ny)->node_pt(l2 + Np * l1) =
            finite_element_pt(Nx * (Ny - 1) + (k - 1) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // We jump to the third dimension
    for (unsigned l3 = 1; l3 < Np; l3++) {
      // The first row of nodes is copied from the element below
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * (Ny - 1) + k * Nx * Ny)
            ->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(Nx * (Ny - 2) + k * Nx * Ny)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      // Second row of  nodes
      // First column of nodes
      for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
        // Determine number of values
        local_node_num = Np * l1 + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = Np * l1 + l2 + l3 * Np * Np;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }

      // Final row of nodes
      // First column of nodes
      // Top left node
      // Determine number of values
      local_node_num = Np * (Np - 1) + l3 * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundaries
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(4, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
        // Allocate memory for the node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        if (Yperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 1);
        }

        // Add the node to the boundary 3
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // Now loop over the rest of the elements in the row, apart from the last
    for (unsigned j = 1; j < (Nx - 1); j++) {
      // Allocate memory for element
      element_num = Nx * (Ny - 1) + j + k * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < Np; l1++) {
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
              ->node_pt(l2 + Np * l1) =
              finite_element_pt(Nx * (Ny - 1) + j + (k - 1) * Nx * Ny)
                  ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
        }
      }

      // Jump in the third dimension

      for (unsigned l3 = 1; l3 < Np; l3++) {
        // The first row is copied from the lower element
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
              ->node_pt(l2 + l3 * Np * Np) =
              finite_element_pt(Nx * (Ny - 2) + j + k * Nx * Ny)
                  ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
        }

        // Second rows
        for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
          // First column is same as neighbouring element
          finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
              ->node_pt(Np * l1 + l3 * Np * Np) =
              finite_element_pt(Nx * (Ny - 1) + (j - 1) + k * Nx * Ny)
                  ->node_pt(Np * l1 + (Np - 1) + l3 * Np * Np);

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < Np; l2++) {
            // Determine number of values
            local_node_num = Np * l1 + l2 + l3 * Np * Np;
            // Allocate memory for the node
            Node_pt[node_count] =
                finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);

            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];

            // Get the fractional position of the node
            finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
            Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (k + s_fraction[2]);

            // No boundaries

            // Increment the node number
            // add_boundary_node(0,Node_pt[node_count]);
            node_count++;
          }
        }

        // Top row
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + j + k * Nx * Ny)
            ->node_pt(Np * (Np - 1) + l3 * Np * Np) =
            finite_element_pt(Nx * (Ny - 1) + (j - 1) + k * Nx * Ny)
                ->node_pt(Np * (Np - 1) + (Np - 1) + l3 * Np * Np);
        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymax;
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          if (Yperiodic) {
            set_periodic_node(Node_pt[node_count], element_num, 1);
          }

          // Add the node to the boundary
          add_boundary_node(3, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }

    } // End of loop over central elements in row

    // FINAL ELEMENT IN ROW (upper right corner)
    //-----------------------------------------

    // Allocate memory for element
    element_num = Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
            ->node_pt(l2 + Np * l1) =
            finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (k - 1) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // We jump to the third dimension

    for (unsigned l3 = 1; l3 < Np; l3++) {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
            ->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(Nx * (Ny - 2) + Nx - 1 + k * Nx * Ny)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
            ->node_pt(Np * l1 + l3 * Np * Np) =
            finite_element_pt(Nx * (Ny - 1) + Nx - 2 + k * Nx * Ny)
                ->node_pt(Np * l1 + (Np - 1) + l3 * Np * Np);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
          // Determine number of values
          local_node_num = Np * l1 + l2 + l3 * Np * Np;
          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];
          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Final node
        // Determine number of values
        local_node_num = Np * l1 + (Np - 1) + l3 * Np * Np;
        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        // if required, make it periodic with the node on the other side
        if (Xperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 0);
        }

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;

      } // End of loop over middle rows

      // Final row
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + k * Nx * Ny)
          ->node_pt(Np * (Np - 1) + l3 * Np * Np) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 2 + k * Nx * Ny)
              ->node_pt(Np * (Np - 1) + (Np - 1) + l3 * Np * Np);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Determine number of values
        local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

        if (Yperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 1);
        }

        // Add the node to the boundary 3
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Final node
      // Determine number of values
      local_node_num = Np * (Np - 1) + (Np - 1) + l3 * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmin + el_length[2] * (k + s_fraction[2]);

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      } else if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundaries 2 and 3
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

  } // End loop of the elements layer

  // END OF THE INTERMEDIATE LAYERS

  // ----------------------------------------------------------------------------------
  //  ****************BEGINNING OF THE LAST
  //  LAYER**************************************
  // ----------------------------------------------------------------------------------

  // FIRST ELEMENT OF THE UPPER LAYER
  //----------------------------------

  element_num = (Nz - 1) * Nx * Ny;
  Element_pt[element_num] = new ELEMENT;

  // The lowest nodes are copied from the lower element
  for (unsigned l1 = 0; l1 < Np; l1++) {
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt((Nz - 1) * Nx * Ny)->node_pt(l2 + Np * l1) =
          finite_element_pt((Nz - 2) * Nx * Ny)
              ->node_pt(l2 + Np * l1 + Np * Np * (Np - 1));
    }
  }

  // Loop over the other node columns in the z direction but the last

  for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
    // Set the corner nodes
    // Determine number of values at this node
    local_node_num = Np * Np * l3;

    // Allocate memory for the node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);

    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];
    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

    // Add the node to boundaries 1 and 4
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Loop over the other nodes in the first row in the x direction
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Determine the number of values at this node
      local_node_num = l2 + Np * Np * l3;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];
      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundary 1
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Loop over the other node columns
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // Determine the number of values
      local_node_num = l1 * Np + Np * Np * l3;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the first node of the row (with boundary 4)
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundary 4
      add_boundary_node(4, Node_pt[node_count]);
      // Increment the node number
      node_count++;

      // Loop over the other nodes in the row
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Set the number of values
        local_node_num = l1 * Np + l2 + Np * Np * l3;

        // Allocate the memory for the node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // No boundary

        // Increment the node number
        node_count++;
      }
    }
  }

  // Make the top nodes

  // Set the corner nodes
  // Determine number of values at this node
  local_node_num = Np * Np * (Np - 1);

  // Allocate memory for the node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);

  // Set the pointer from the element
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmin;
  Node_pt[node_count]->x(1) = Ymin;
  Node_pt[node_count]->x(2) = Zmax;

  // Add the node to boundaries 1, 4 and 5
  add_boundary_node(1, Node_pt[node_count]);
  add_boundary_node(4, Node_pt[node_count]);
  add_boundary_node(5, Node_pt[node_count]);

  // Increment the node number
  node_count++;

  // Loop over the other nodes in the first row in the x direction
  for (unsigned l2 = 1; l2 < Np; l2++) {
    // Determine the number of values at this node
    local_node_num = l2 + Np * Np * (Np - 1);

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundaries 1 and 5
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);
    // Increment the node number
    node_count++;
  }

  // Loop over the other node columns
  for (unsigned l1 = 1; l1 < Np; l1++) {
    // Determine the number of values
    local_node_num = l1 * Np + Np * Np * (Np - 1);

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the first node of the row (with boundary 4)
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundary 4
    add_boundary_node(4, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Loop over the other nodes in the row
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Set the number of values
      local_node_num = l1 * Np + l2 + Np * Np * (Np - 1);

      // Allocate the memory for the node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundary 5
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }
  }

  //----------------------------------------------------------------------

  // CENTRE OF FIRST ROW OF ELEMENTS OF THE UPPER LAYER
  //--------------------------------------------------

  // Now loop over the first row of elements, apart from final element
  for (unsigned j = 1; j < (Nx - 1); j++) {
    // Allocate memory for new element
    element_num = j + (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(j + (Nz - 1) * Nx * Ny)->node_pt(l2 + Np * l1) =
            finite_element_pt(j + (Nz - 2) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // Loop over the other node columns in the z direction but the last
    for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(j + (Nz - 1) * Nx * Ny)->node_pt(l3 * Np * Np) =
          finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
              ->node_pt(l3 * Np * Np + (Np - 1));

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = l2 + l3 * Np * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin;
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 1
        add_boundary_node(1, Node_pt[node_count]);
        // Increment the node number
        node_count++;
      }

      // Do the rest of the nodes
      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column of nodes is same as neighbouring element
        finite_element_pt(j + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * Np + l3 * Np * Np) =
            finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
                ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = l1 * Np + l2 + l3 * Np * Np;

          // Allocate memory for the nodes
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
          Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }
    }

    // Top nodes

    // First node is same as neighbouring element
    finite_element_pt(j + (Nz - 1) * Nx * Ny)->node_pt((Np - 1) * Np * Np) =
        finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt((Np - 1) * Np * Np + (Np - 1));

    // New nodes for other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Determine number of values
      local_node_num = l2 + (Np - 1) * Np * Np;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 1 and 5
      add_boundary_node(1, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Do the rest of the nodes
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(j + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * Np + (Np - 1) * Np * Np) =
          finite_element_pt(j - 1 + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * Np + (Np - 1) + (Np - 1) * Np * Np);

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = l1 * Np + l2 + (Np - 1) * Np * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) = Zmax;

        // Add to the boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }
  } // End loop of first row of elements

  // MY FINAL ELEMENT IN THE FIRST ROW OF THE UPPER LAYER (right corner)
  //--------------------------------------------------------------------

  // Allocate memory for new element
  element_num = Nx - 1 + (Nz - 1) * Nx * Ny;
  Element_pt[element_num] = new ELEMENT;

  // The lowest nodes are copied from the lower element
  for (unsigned l1 = 0; l1 < Np; l1++) {
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)->node_pt(l2 + Np * l1) =
          finite_element_pt(Nx - 1 + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
    }
  }

  // Loop over the other node columns in the z direction but the last
  for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
    // First column of nodes is same as neighbouring element
    finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)->node_pt(l3 * Np * Np) =
        finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(l3 * Np * Np + (Np - 1));

    // New nodes for other rows (y=0)
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Determine number of values
      local_node_num = l2 + l3 * Np * Np;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin;
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundary 1
      add_boundary_node(1, Node_pt[node_count]);
      // Increment the node number
      node_count++;
    }

    // Last node of the row (y=0)

    // Determine number of values
    local_node_num = (Np - 1) + l3 * Np * Np;

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    }

    // Add the node to the boundary 1 and 2
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(2, Node_pt[node_count]);
    // Increment the node number
    node_count++;

    // Do the rest of the nodes
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column of nodes is same as neighbouring element
      finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * Np + l3 * Np * Np) =
          finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Determine number of values
        local_node_num = l1 * Np + l2 + l3 * Np * Np;

        // Allocate memory for the nodes
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // No boundaries

        // Increment the node number
        node_count++;
      }

      // Last nodes (at the surface x = Lx)
      // Determine number of values
      local_node_num = l1 * Np + (Np - 1) + l3 * Np * Np;
      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundary 2
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }
  }

  // We make the top nodes
  // First node is same as neighbouring element
  finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)->node_pt((Np - 1) * Np * Np) =
      finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
          ->node_pt((Np - 1) * Np * Np + (Np - 1));

  // New nodes for other rows (y=0)
  for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
    // Determine number of values
    local_node_num = l2 + (Np - 1) * Np * Np;

    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
    Node_pt[node_count]->x(1) = Ymin;
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundaries 1 and 5
    add_boundary_node(1, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // Last node of the row (y=0)

  // Determine number of values
  local_node_num = (Np - 1) + (Np - 1) * Np * Np;

  // Allocate memory for the nodes
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer from the element
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmax;
  Node_pt[node_count]->x(1) = Ymin;
  Node_pt[node_count]->x(2) = Zmax;

  // if required, make it periodic with the node on the other side
  if (Xperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 0);
  }

  // Add the node to the boundary 1 and 2
  add_boundary_node(1, Node_pt[node_count]);
  add_boundary_node(2, Node_pt[node_count]);
  add_boundary_node(5, Node_pt[node_count]);
  // Increment the node number
  node_count++;

  // Do the rest of the nodes
  for (unsigned l1 = 1; l1 < Np; l1++) {
    // First column of nodes is same as neighbouring element
    finite_element_pt(Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(l1 * Np + (Np - 1) * Np * Np) =
        finite_element_pt(Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * Np + (Np - 1) + (Np - 1) * Np * Np);

    // New nodes for other columns
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Determine number of values
      local_node_num = l1 * Np + l2 + (Np - 1) * Np * Np;

      // Allocate memory for the nodes
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
      Node_pt[node_count]->x(2) = Zmax;

      // Add node to the boundary 5
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Last nodes (at the surface x = Lx)
    // Determine number of values
    local_node_num = l1 * Np + (Np - 1) + (Np - 1) * Np * Np;
    // Allocate memory for the nodes
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * s_fraction[1];
    Node_pt[node_count]->x(2) = Zmax;

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    }

    // Add the node to the boundaries 2 and 5
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // ALL CENTRAL ELEMENT ROWS OF THE TOP  LAYER
  //------------------------------------------

  // Loop over remaining element rows
  for (unsigned i = 1; i < (Ny - 1); i++) {
    // Set the first element in the row

    // Allocate memory for element (x =0)
    element_num = Nx * i + Nx * Ny * (Nz - 1);
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * i + (Nz - 1) * Nx * Ny)->node_pt(l2 + Np * l1) =
            finite_element_pt(Nx * i + (Nz - 2) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // We extend now this element to the third dimension

    for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
      // The first row of nodes is copied from the element below (at z=0)
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * i + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(Nx * (i - 1) + (Nz - 1) * Nx * Ny)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      // Other rows are new nodes (first the nodes for which x=0)
      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column of nodes

        // Determine number of values
        local_node_num = l1 * Np + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // Add the node to the boundary 4
        add_boundary_node(4, Node_pt[node_count]);

        // Increment the node number
        node_count++;

        // Now do the other columns (we extend to the rest of the nodes)
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = l1 * Np + l2 + Np * Np * l3;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer from the element
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }
      }
    }

    // Now the top nodes

    // The first row of nodes is copied from the element below (at z=0)
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * i + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + (Np - 1) * Np * Np) =
          finite_element_pt(Nx * (i - 1) + (Nz - 1) * Nx * Ny)
              ->node_pt((Np - 1) * Np + l2 + (Np - 1) * Np * Np);
    }

    // Other rows are new nodes (first the nodes for which x=0)
    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column of nodes

      // Determine number of values
      local_node_num = l1 * Np + (Np - 1) * Np * Np;

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundaries 4 and 5
      add_boundary_node(4, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns (we extend to the rest of the nodes)
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = l1 * Np + l2 + Np * Np * (Np - 1);

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // Now loop over the rest of the elements in the row, apart from the last
    for (unsigned j = 1; j < (Nx - 1); j++) {
      // Allocate memory for new element
      element_num = Nx * i + j + (Nz - 1) * Nx * Ny;
      Element_pt[element_num] = new ELEMENT;

      // The lowest nodes are copied from the lower element
      for (unsigned l1 = 0; l1 < Np; l1++) {
        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
              ->node_pt(l2 + Np * l1) =
              finite_element_pt(Nx * i + j + (Nz - 2) * Nx * Ny)
                  ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
        }
      }

      // We extend to the third dimension but the last layer of nodes

      for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
        // The first row is copied from the lower element

        for (unsigned l2 = 0; l2 < Np; l2++) {
          finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
              ->node_pt(l2 + l3 * Np * Np) =
              finite_element_pt(Nx * (i - 1) + j + (Nz - 1) * Nx * Ny)
                  ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
        }

        for (unsigned l1 = 1; l1 < Np; l1++) {
          // First column is same as neighbouring element
          finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * Np + l3 * Np * Np) =
              finite_element_pt(Nx * i + (j - 1) + (Nz - 1) * Nx * Ny)
                  ->node_pt(l1 * Np + l3 * Np * Np + (Np - 1));

          // New nodes for other columns
          for (unsigned l2 = 1; l2 < Np; l2++) {
            // Determine number of values
            local_node_num = l1 * Np + l2 + l3 * Np * Np;

            // Allocate memory for the nodes
            Node_pt[node_count] =
                finite_element_pt(element_num)
                    ->construct_node(local_node_num, time_stepper_pt);
            // Set the pointer
            finite_element_pt(element_num)->node_pt(local_node_num) =
                Node_pt[node_count];
            // Get the fractional position of the node
            finite_element_pt(element_num)
                ->local_fraction_of_node(local_node_num, s_fraction);

            // Set the position of the node
            Node_pt[node_count]->x(0) =
                Xmin + el_length[0] * (j + s_fraction[0]);
            Node_pt[node_count]->x(1) =
                Ymin + el_length[1] * (i + s_fraction[1]);
            Node_pt[node_count]->x(2) =
                Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);
            // No boundaries

            // Increment the node number
            node_count++;
          }
        }
      }

      // Now the top nodes

      // The first row is copied from the lower element

      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + (Np - 1) * Np * Np) =
            finite_element_pt(Nx * (i - 1) + j + (Nz - 1) * Nx * Ny)
                ->node_pt((Np - 1) * Np + l2 + (Np - 1) * Np * Np);
      }

      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column is same as neighbouring element
        finite_element_pt(Nx * i + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * Np + (Np - 1) * Np * Np) =
            finite_element_pt(Nx * i + (j - 1) + (Nz - 1) * Nx * Ny)
                ->node_pt(l1 * Np + (Np - 1) * Np * Np + (Np - 1));

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = l1 * Np + l2 + (Np - 1) * Np * Np;

          // Allocate memory for the nodes
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_boundary_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) = Zmax;

          // Add to boundary 5
          add_boundary_node(5, Node_pt[node_count]);

          // Increment the node number
          node_count++;
        }
      }

    } // End of loop over elements in row

    // Do final element in row

    // Allocate memory for element
    element_num = Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + Np * l1) =
            finite_element_pt(Nx * i + Nx - 1 + (Nz - 2) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // We go to the third dimension but we do not reach the top
    for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(Nx * (i - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      for (unsigned l1 = 1; l1 < Np; l1++) {
        // First column is same as neighbouring element
        finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt(l1 * Np + l3 * Np * Np) =
            finite_element_pt(Nx * i + Nx - 2 + (Nz - 1) * Nx * Ny)
                ->node_pt(l1 * Np + (Np - 1) + l3 * Np * Np);

        // Middle nodes
        for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
          // Determine number of values
          local_node_num = l1 * Np + l2 + l3 * Np * Np;

          // Allocate memory for node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);
          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) =
              Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
          Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
          Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          node_count++;
        }

        // Final node

        // Determine number of values
        local_node_num = l1 * Np + (Np - 1) + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmax;
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // if required, make it periodic with the node on the other side
        if (Xperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 0);
        }

        // Add the node to the boundary 2
        add_boundary_node(2, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      } // End of loop over rows of nodes in the element

    } // End of the 3dimension loop

    // Now the top nodes

    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + (Np - 1) * Np * Np) =
          finite_element_pt(Nx * (i - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
              ->node_pt((Np - 1) * Np + l2 + (Np - 1) * Np * Np);
    }

    for (unsigned l1 = 1; l1 < Np; l1++) {
      // First column is same as neighbouring element
      finite_element_pt(Nx * i + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l1 * Np + (Np - 1) * Np * Np) =
          finite_element_pt(Nx * i + Nx - 2 + (Nz - 1) * Nx * Ny)
              ->node_pt(l1 * Np + (Np - 1) + (Np - 1) * Np * Np);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Determine number of values
        local_node_num = l1 * Np + l2 + (Np - 1) * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add to boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }

      // Final node

      // Determine number of values
      local_node_num = l1 * Np + (Np - 1) + (Np - 1) * Np * Np;

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) = Ymin + el_length[1] * (i + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmax;

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundaries 2 and 5
      add_boundary_node(2, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    } // End of loop over rows of nodes in the element

  } // End of loop over rows of elements

  // FINAL ELEMENT ROW (IN TOP  LAYER)
  //===========================================

  // FIRST ELEMENT IN UPPER ROW (upper left corner)
  //----------------------------------------------

  // Allocate memory for element
  element_num = Nx * (Ny - 1) + (Nz - 1) * Nx * Ny;
  Element_pt[element_num] = new ELEMENT;

  // The lowest nodes are copied from the lower element
  for (unsigned l1 = 0; l1 < Np; l1++) {
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * (Ny - 1) + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + Np * l1) =
          finite_element_pt(Nx * (Ny - 1) + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
    }
  }

  // We jump to the third dimension but the last layer of nodes
  for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
    // The first row of nodes is copied from the element below
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * (Ny - 1) + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + l3 * Np * Np) =
          finite_element_pt(Nx * (Ny - 2) + (Nz - 1) * Nx * Ny)
              ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
    }

    // Second row of  nodes
    // First column of nodes
    for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
      // Determine number of values
      local_node_num = Np * l1 + l3 * Np * Np;

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin;
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // Add the node to the boundary 4
      add_boundary_node(4, Node_pt[node_count]);

      // Increment the node number
      node_count++;

      // Now do the other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = Np * l1 + l2 + l3 * Np * Np;

        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer from the element
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // No boundaries

        // Increment the node number
        node_count++;
      }
    }

    // Final row of nodes
    // First column of nodes
    // Top left node
    // Determine number of values
    local_node_num = Np * (Np - 1) + l3 * Np * Np;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

    if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(4, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Determine number of values
      local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
      // Allocate memory for the node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundary 3
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }
  }
  // Now the top nodes

  // The first row of nodes is copied from the element below
  for (unsigned l2 = 0; l2 < Np; l2++) {
    finite_element_pt(Nx * (Ny - 1) + (Nz - 1) * Nx * Ny)
        ->node_pt(l2 + (Np - 1) * Np * Np) =
        finite_element_pt(Nx * (Ny - 2) + (Nz - 1) * Nx * Ny)
            ->node_pt((Np - 1) * Np + l2 + (Np - 1) * Np * Np);
  }

  // Second row of  nodes
  // First column of nodes
  for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
    // Determine number of values
    local_node_num = Np * l1 + (Np - 1) * Np * Np;

    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
    Node_pt[node_count]->x(2) = Zmax;

    // Add the node to the boundaries 4 and 5
    add_boundary_node(4, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;

    // Now do the other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Determine number of values
      local_node_num = Np * l1 + l2 + (Np - 1) * Np * Np;

      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer from the element
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];
      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmax;

      // Add the node to the boundary 5
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }
  }

  // Final row of nodes
  // First column of nodes
  // Top left node
  // Determine number of values
  local_node_num = Np * (Np - 1) + (Np - 1) * Np * Np;
  // Allocate memory for node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer from the element
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];
  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmin;
  Node_pt[node_count]->x(1) = Ymax;
  Node_pt[node_count]->x(2) = Zmax;

  if (Yperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 1);
  }

  // Add the node to the boundaries
  add_boundary_node(3, Node_pt[node_count]);
  add_boundary_node(4, Node_pt[node_count]);
  add_boundary_node(5, Node_pt[node_count]);

  // Increment the node number
  node_count++;

  // Now do the other columns
  for (unsigned l2 = 1; l2 < Np; l2++) {
    // Determine number of values
    local_node_num = Np * (Np - 1) + l2 + (Np - 1) * Np * Np;
    // Allocate memory for the node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer from the element
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * s_fraction[0];
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmax;

    if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries 3 and 5
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // Now loop over the rest of the elements in the row, apart from the last
  for (unsigned j = 1; j < (Nx - 1); j++) {
    // Allocate memory for element
    element_num = Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny;
    Element_pt[element_num] = new ELEMENT;

    // The lowest nodes are copied from the lower element
    for (unsigned l1 = 0; l1 < Np; l1++) {
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + Np * l1) =
            finite_element_pt(Nx * (Ny - 1) + j + (Nz - 2) * Nx * Ny)
                ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
      }
    }

    // Jump in the third dimension but the top nodes

    for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
      // The first row is copied from the lower element
      for (unsigned l2 = 0; l2 < Np; l2++) {
        finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
            ->node_pt(l2 + l3 * Np * Np) =
            finite_element_pt(Nx * (Ny - 2) + j + (Nz - 1) * Nx * Ny)
                ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
      }

      // Second rows
      for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
        // First column is same as neighbouring element
        finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
            ->node_pt(Np * l1 + l3 * Np * Np) =
            finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
                ->node_pt(Np * l1 + (Np - 1) + l3 * Np * Np);

        // New nodes for other columns
        for (unsigned l2 = 1; l2 < Np; l2++) {
          // Determine number of values
          local_node_num = Np * l1 + l2 + l3 * Np * Np;
          // Allocate memory for the node
          Node_pt[node_count] =
              finite_element_pt(element_num)
                  ->construct_node(local_node_num, time_stepper_pt);

          // Set the pointer
          finite_element_pt(element_num)->node_pt(local_node_num) =
              Node_pt[node_count];

          // Get the fractional position of the node
          finite_element_pt(element_num)
              ->local_fraction_of_node(local_node_num, s_fraction);

          // Set the position of the node
          Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
          Node_pt[node_count]->x(1) =
              Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
          Node_pt[node_count]->x(2) =
              Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

          // No boundaries

          // Increment the node number
          // add_boundary_node(0,Node_pt[node_count]);
          node_count++;
        }
      }

      // Top row
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
          ->node_pt(Np * (Np - 1) + l3 * Np * Np) =
          finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
              ->node_pt(Np * (Np - 1) + (Np - 1) + l3 * Np * Np);
      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) = Ymax;
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        if (Yperiodic) {
          set_periodic_node(Node_pt[node_count], element_num, 1);
        }

        // Add the node to the boundary
        add_boundary_node(3, Node_pt[node_count]);

        // Increment the node number
        node_count++;
      }
    }

    // Now the top nodes

    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + (Np - 1) * Np * Np) =
          finite_element_pt(Nx * (Ny - 2) + j + (Nz - 1) * Nx * Ny)
              ->node_pt((Np - 1) * Np + l2 + (Np - 1) * Np * Np);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
          ->node_pt(Np * l1 + (Np - 1) * Np * Np) =
          finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
              ->node_pt(Np * l1 + (Np - 1) + (Np - 1) * Np * Np);

      // New nodes for other columns
      for (unsigned l2 = 1; l2 < Np; l2++) {
        // Determine number of values
        local_node_num = Np * l1 + l2 + (Np - 1) * Np * Np;
        // Allocate memory for the node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_boundary_node(local_node_num, time_stepper_pt);

        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) = Zmax;

        // Add the node to the boundary 5
        add_boundary_node(5, Node_pt[node_count]);

        // Increment the node number add_boundary_node(0,Node_pt[node_count]);
        node_count++;
      }
    }

    // Top row
    // First column is same as neighbouring element
    finite_element_pt(Nx * (Ny - 1) + j + (Nz - 1) * Nx * Ny)
        ->node_pt(Np * (Np - 1) + (Np - 1) * Np * Np) =
        finite_element_pt(Nx * (Ny - 1) + (j - 1) + (Nz - 1) * Nx * Ny)
            ->node_pt(Np * (Np - 1) + (Np - 1) + (Np - 1) * Np * Np);
    // New nodes for other columns
    for (unsigned l2 = 1; l2 < Np; l2++) {
      // Determine number of values
      local_node_num = Np * (Np - 1) + l2 + (Np - 1) * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmin + el_length[0] * (j + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) = Zmax;

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundaries 3 and 5
      add_boundary_node(3, Node_pt[node_count]);
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

  } // End of loop over central elements in row

  // LAST ELEMENT (Really!!)
  //-----------------------------------------

  // Allocate memory for element
  element_num = Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny;
  Element_pt[element_num] = new ELEMENT;

  // The lowest nodes are copied from the lower element
  for (unsigned l1 = 0; l1 < Np; l1++) {
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + Np * l1) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 2) * Nx * Ny)
              ->node_pt(l2 + Np * l1 + (Np - 1) * Np * Np);
    }
  }

  // We jump to the third dimension but the top nodes

  for (unsigned l3 = 1; l3 < (Np - 1); l3++) {
    // The first row is copied from the lower element
    for (unsigned l2 = 0; l2 < Np; l2++) {
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(l2 + l3 * Np * Np) =
          finite_element_pt(Nx * (Ny - 2) + Nx - 1 + (Nz - 1) * Nx * Ny)
              ->node_pt((Np - 1) * Np + l2 + l3 * Np * Np);
    }

    // Second rows
    for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
      // First column is same as neighbouring element
      finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
          ->node_pt(Np * l1 + l3 * Np * Np) =
          finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
              ->node_pt(Np * l1 + (Np - 1) + l3 * Np * Np);

      // Middle nodes
      for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
        // Determine number of values
        local_node_num = Np * l1 + l2 + l3 * Np * Np;
        // Allocate memory for node
        Node_pt[node_count] =
            finite_element_pt(element_num)
                ->construct_node(local_node_num, time_stepper_pt);
        // Set the pointer
        finite_element_pt(element_num)->node_pt(local_node_num) =
            Node_pt[node_count];

        // Get the fractional position of the node
        finite_element_pt(element_num)
            ->local_fraction_of_node(local_node_num, s_fraction);

        // Set the position of the node
        Node_pt[node_count]->x(0) =
            Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
        Node_pt[node_count]->x(1) =
            Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
        Node_pt[node_count]->x(2) =
            Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

        // No boundaries

        // Increment the node number
        node_count++;
      }

      // Final node
      // Determine number of values
      local_node_num = Np * l1 + (Np - 1) + l3 * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) = Xmax;
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      // if required, make it periodic with the node on the other side
      if (Xperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 0);
      }

      // Add the node to the boundary 2
      add_boundary_node(2, Node_pt[node_count]);

      // Increment the node number
      node_count++;

    } // End of loop over middle rows

    // Final row
    // First column is same as neighbouring element
    finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(Np * (Np - 1) + l3 * Np * Np) =
        finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(Np * (Np - 1) + (Np - 1) + l3 * Np * Np);

    // Middle nodes
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Determine number of values
      local_node_num = Np * (Np - 1) + l2 + l3 * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) = Ymax;
      Node_pt[node_count]->x(2) =
          Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

      if (Yperiodic) {
        set_periodic_node(Node_pt[node_count], element_num, 1);
      }

      // Add the node to the boundary 3
      add_boundary_node(3, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Final node
    // Determine number of values
    local_node_num = Np * (Np - 1) + (Np - 1) + l3 * Np * Np;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    // In fluid 2
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmin + el_length[2] * (Nz - 1 + s_fraction[2]);

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    } else if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundaries 2 and 3
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(3, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // Now the top nodes

  // The first row is copied from the lower element
  for (unsigned l2 = 0; l2 < Np; l2++) {
    finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(l2 + (Np - 1) * Np * Np) =
        finite_element_pt(Nx * (Ny - 2) + Nx - 1 + (Nz - 1) * Nx * Ny)
            ->node_pt((Np - 1) * Np + l2 + (Np - 1) * Np * Np);
  }

  // Second rows
  for (unsigned l1 = 1; l1 < (Np - 1); l1++) {
    // First column is same as neighbouring element
    finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
        ->node_pt(Np * l1 + (Np - 1) * Np * Np) =
        finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
            ->node_pt(Np * l1 + (Np - 1) + (Np - 1) * Np * Np);

    // Middle nodes
    for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
      // Determine number of values
      local_node_num = Np * l1 + l2 + (Np - 1) * Np * Np;
      // Allocate memory for node
      Node_pt[node_count] =
          finite_element_pt(element_num)
              ->construct_boundary_node(local_node_num, time_stepper_pt);
      // Set the pointer
      finite_element_pt(element_num)->node_pt(local_node_num) =
          Node_pt[node_count];

      // Get the fractional position of the node
      finite_element_pt(element_num)
          ->local_fraction_of_node(local_node_num, s_fraction);

      // Set the position of the node
      Node_pt[node_count]->x(0) =
          Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
      Node_pt[node_count]->x(1) =
          Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
      Node_pt[node_count]->x(2) = Zmax;

      // Add to boundary 5
      add_boundary_node(5, Node_pt[node_count]);

      // Increment the node number
      node_count++;
    }

    // Final node
    // Determine number of values
    local_node_num = Np * l1 + (Np - 1) + (Np - 1) * Np * Np;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmax;
    Node_pt[node_count]->x(1) = Ymin + el_length[1] * (Ny - 1 + s_fraction[1]);
    Node_pt[node_count]->x(2) = Zmax;

    // if required, make it periodic with the node on the other side
    if (Xperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 0);
    }

    // Add the node to the boundaries 2 and 5
    add_boundary_node(2, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;

  } // End of loop over middle rows

  // Final row
  // First node is same as neighbouring element
  finite_element_pt(Nx * (Ny - 1) + Nx - 1 + (Nz - 1) * Nx * Ny)
      ->node_pt(Np * (Np - 1) + (Np - 1) * Np * Np) =
      finite_element_pt(Nx * (Ny - 1) + Nx - 2 + (Nz - 1) * Nx * Ny)
          ->node_pt(Np * (Np - 1) + (Np - 1) + (Np - 1) * Np * Np);

  // Middle nodes
  for (unsigned l2 = 1; l2 < (Np - 1); l2++) {
    // Determine number of values
    local_node_num = Np * (Np - 1) + l2 + (Np - 1) * Np * Np;
    // Allocate memory for node
    Node_pt[node_count] =
        finite_element_pt(element_num)
            ->construct_boundary_node(local_node_num, time_stepper_pt);
    // Set the pointer
    finite_element_pt(element_num)->node_pt(local_node_num) =
        Node_pt[node_count];

    // Get the fractional position of the node
    finite_element_pt(element_num)
        ->local_fraction_of_node(local_node_num, s_fraction);

    // Set the position of the node
    Node_pt[node_count]->x(0) = Xmin + el_length[0] * (Nx - 1 + s_fraction[0]);
    Node_pt[node_count]->x(1) = Ymax;
    Node_pt[node_count]->x(2) = Zmax;

    if (Yperiodic) {
      set_periodic_node(Node_pt[node_count], element_num, 1);
    }

    // Add the node to the boundary 3
    add_boundary_node(3, Node_pt[node_count]);
    add_boundary_node(5, Node_pt[node_count]);

    // Increment the node number
    node_count++;
  }

  // Final node (really!!)
  // Determine number of values
  local_node_num = Np * (Np - 1) + (Np - 1) + (Np - 1) * Np * Np;
  // Allocate memory for node
  Node_pt[node_count] =
      finite_element_pt(element_num)
          ->construct_boundary_node(local_node_num, time_stepper_pt);
  // Set the pointer
  finite_element_pt(element_num)->node_pt(local_node_num) = Node_pt[node_count];

  // Get the fractional position of the node
  finite_element_pt(element_num)
      ->local_fraction_of_node(local_node_num, s_fraction);

  // Set the position of the node
  Node_pt[node_count]->x(0) = Xmax;
  Node_pt[node_count]->x(1) = Ymax;
  Node_pt[node_count]->x(2) = Zmax;

  // if required, make it periodic with the node on the other side
  if (Xperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 0);
  } else if (Yperiodic) {
    set_periodic_node(Node_pt[node_count], element_num, 1);
  }

  // Add the node to the boundaries 2, 3 and 5
  add_boundary_node(2, Node_pt[node_count]);
  add_boundary_node(3, Node_pt[node_count]);
  add_boundary_node(5, Node_pt[node_count]);

  // Increment the node number
  node_count++;

  // Setup lookup scheme that establishes which elements are located
  // on the mesh boundaries
  setup_boundary_element_info();
}

#endif // CUBIC_BRICK_MESH_H
