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
  // Build the single layer mesh
  build_single_layer_mesh(time_stepper_pt);
}


template <class ELEMENT>
SingleLayerSpineMesh3D<ELEMENT>::SingleLayerSpineMesh3D(
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &lx, const double &ly, const double &h,
    const bool &periodic_in_x, const bool &periodic_in_y,
    TimeStepper *time_stepper_pt)
    : CubicBrickMesh<ELEMENT>(nx, ny, nz, lx, ly, h, periodic_in_x, periodic_in_y,
                              time_stepper_pt) {
  // Build the single layer mesh
  build_single_layer_mesh(time_stepper_pt);
}


template <class ELEMENT>
void SingleLayerSpineMesh3D<ELEMENT>::build_single_layer_mesh(
    TimeStepper *time_stepper_pt) {
  // mesh can only be built with 3D Qelements
  MeshChecker::assert_geometric_element<QElementGeometricBase, ELEMENT>(3);

  // build the underlying brick mesh
  CubicBrickMesh<ELEMENT>::build_mesh(time_stepper_pt);

  // TODO: deal with the spines
}

#endif // SINGLELAYERSPINEMESH3D_H
