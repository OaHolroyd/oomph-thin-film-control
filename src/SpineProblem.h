//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef SPINEINCLINEDPLANEPROBLEM_H
#define SPINEINCLINEDPLANEPROBLEM_H

#include "fluid_interface.h"
#include "navier_stokes.h"

#include "Problem.h"
#include "single_layer_cubic_spine_mesh.h"

// ======================================================================
//   Create a spine mesh for the problem
// ======================================================================
template <class ELEMENT>
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
  SpineFilmMesh(const unsigned &nx, const unsigned &ny, const unsigned &nz,
                const double &lx, const double &ly, const double &lz,
                TimeStepper *time_stepper_pt)
      : SingleLayerCubicSpineMesh<ELEMENT>(nx, ny, nz, lx, ly, lz, true, true,
                                           time_stepper_pt) {
  } // end of constructor
};

// ============================================================================
//  Specific class for controlled film problem using spines
// ============================================================================
template <class ELEMENT, class TIMESTEPPER>
class SpineControlledFilmProblem
    : public ControlledFilmProblem<ELEMENT,
                                   SpineSurfaceFluidInterfaceElement<ELEMENT>> {
public:
  /**
   * Constructor for the spine controlled film problem
   *
   * @param nx Number of elements in the x direction
   * @param ny Number of elements in the y direction
   * @param nz Number of elements in the z direction
   * @param nx_control Number of elements in the x direction of the control
   * system
   * @param ny_control Number of elements in the y direction of the control
   * system
   * @param m_control the number of actuators
   * @param p_control the number of observers
   */
  SpineControlledFilmProblem(const unsigned &nx, const unsigned &ny,
                             const unsigned &nz, const int &nx_control,
                             const int &ny_control, const int &m_control,
                             const int &p_control, const bool &use_mumps)
      : ControlledFilmProblem<ELEMENT,
                              SpineSurfaceFluidInterfaceElement<ELEMENT>>(
            nx_control, ny_control, m_control, p_control, use_mumps) {
    using namespace Global_Variables; // to access the length

    // create our one and only timestepper, with adaptive timestepping
    this->add_time_stepper_pt(new TIMESTEPPER);

    // create the bulk mesh
    this->Bulk_mesh_pt = new SpineFilmMesh<ELEMENT>(nx, ny, nz, Lx, Ly, 1.0,
                                                    this->time_stepper_pt());
    this->nx = nx;
    this->ny = ny;

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

  /// override the set_hqf function to use information from the Spine Mesh
  void set_hqf(int use_control);

private:
  int nx, ny;
};

template <class ELEMENT, class INTERFACE_ELEMENT>
void SpineControlledFilmProblem<ELEMENT, INTERFACE_ELEMENT>::set_hqf(
    int use_control) {
  using namespace Global_Variables;

  // grid spacings
  double dx = Lx / this->nx_control;
  double dy = Ly / this->ny_control;

  // get values of h in a nx_control x ny_control grid
  for (int i = 0; i < this->ny_control; i++) {
    double yi = (dy * (static_cast<double>(i) + 0.5)); // y coordinate
    double pi = yi / Ly; // normalised y coordinate
    int ei = pi * ny;
    for (int j = 0; j < this->nx_control; j++) {
      double xj = (dx * (static_cast<double>(j) + 0.5)); // x coordinate
      double pj = xj / Lx; // normalised x coordinate
      int ej = pj * nx;
      int k = j + i * this->nx_control; // linear index into the arrays

      // find the surface element containing (xj, yi)
      int e = ej + ei * nx;
      FaceElement *el =
          dynamic_cast<FaceElement *>(this->Surface_mesh_pt->element_pt(e));

      int np = el->nnode_1d();

      // find the nodes in the element that enclose the point
      SpineNode *n0 = dynamic_cast<SpineNode *>(el->node_pt(0));
      SpineNode *n1 = dynamic_cast<SpineNode *>(el->node_pt(np - 1));
      SpineNode *n2 = dynamic_cast<SpineNode *>(el->node_pt((np - 1) * np));
      SpineNode *n3 = dynamic_cast<SpineNode *>(el->node_pt(np * np - 1));

      // find enclosing nodes
      int nnodes = el->nnode();
      for (int n = 0; n < nnodes; n++) {
        SpineNode *nn = dynamic_cast<SpineNode *>(el->node_pt(n));
        double x = nn->x(0);
        double y = nn->x(1);

        if ((x <= xj) && (y <= yi)) {
          // candidate for lower left node
          if (n0->x(0) <= x && n0->x(1) <= y) {
            n0 = nn;
          }
        }

        if ((x >= xj) && (y <= yi)) {
          // candidate for lower right node
          if (n1->x(0) >= x && n1->x(1) <= y) {
            n1 = nn;
          }
        }

        if ((x <= xj) && (y >= yi)) {
          // candidate for upper left node
          if (n2->x(0) <= x && n2->x(1) >= y) {
            n2 = nn;
          }
        }

        if ((x >= xj) && (y >= yi)) {
          // candidate for upper right node
          if (n3->x(0) >= x && n3->x(1) >= y) {
            n3 = nn;
          }
        }
      }

      assert(n0->x(0) <= xj && n0->x(1) <= yi);
      assert(n1->x(0) >= xj && n1->x(1) <= yi);
      assert(n2->x(0) <= xj && n2->x(1) >= yi);
      assert(n3->x(0) >= xj && n3->x(1) >= yi);

      // convert coordinates to local coordinates
      double x0 = n0->x(0);
      double y0 = n0->x(1);
      double x1 = n3->x(0);
      double y1 = n3->x(1);

      double x = 0.5;
      if (x0 != x1) {
        x = (xj - x0) / (x1 - x0);
      }

      double y = 0.5;
      if (y0 != y1) {
        y = (yi - y0) / (y1 - y0);
      }

      // interpolate the height from the four spines to the point
      double h0 = n0->spine_pt()->height();
      double h1 = n1->spine_pt()->height();
      double h2 = n2->spine_pt()->height();
      double h3 = n3->spine_pt()->height();
      this->h[k] = h0 * (1 - x) * (1 - y) + h1 * x * (1 - y) +
                   h2 * (1 - x) * y + h3 * x * y;

      // this->f[k] = (use_control > 0) ? control(xj, yi) : 0.0;
      this->f[k] = (use_control > 0) ? (1.0 - this->h[k]) : 0.0;


      // TODO: qx, and qy
      this->qx[k] = 0.0;
      this->qy[k] = 0.0;
    }
  }
}

#endif // SPINEINCLINEDPLANEPROBLEM_H
