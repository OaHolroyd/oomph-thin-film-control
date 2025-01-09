//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef GLOBALPHYSICALVARIABLES_H
#define GLOBALPHYSICALVARIABLES_H

//Standard C++ library includes
#include <iostream>
#include <cmath>

//Finite-Element library routines
#include "generic.h"
#include "solid.h"

using namespace std;

using namespace oomph;

//The global physical variables
namespace Global_Physical_Variables {
  /// Reynolds number, based on the average velocity within the fluid film
  double Re = 0.0;

  /// The product of Reynolds number and inverse Froude number
  /// is set to two in this problem, which gives the free surface velocity
  /// to be sin(alpha). [Set to three in order to get the same scale as
  /// used by Yih, Benjamin, etc]
  double ReInvFr = 2.0;

  /// Angle of incline of the slope (45 degrees)
  double Alpha = 1.0 * atan(1.0);

  /// The Vector direction of gravity, set in main()
  Vector<double> G(2, 0.0);

  /// The Capillary number
  double Ca = 1.0;

  /// Set the wavenumber
  double K = 0.1;

  /// Set the number of waves desired in the domain
  double N_wave = 3;

  /// The length of the domain to fit the desired number of waves
  double Length = 2 * N_wave * 4.0 * atan(1.0) / K;

  /// Direction of the wall normal vector (at the inlet)
  Vector<double> Wall_normal;

  /// Function that specifies the wall unit normal at the inlet
  void wall_unit_normal_inlet_fct(const Vector<double> &x, Vector<double> &normal) {
    normal = Wall_normal;
  }

  /// Function that specified the wall unit normal at the outlet
  void wall_unit_normal_outlet_fct(const Vector<double> &x, Vector<double> &normal) {
    //Set the normal
    normal = Wall_normal;
    //and flip the sign
    unsigned n_dim = normal.size();
    for (unsigned i = 0; i < n_dim; ++i) { normal[i] *= -1.0; }
  }

  /// The contact angle that is imposed at the inlet (pi)
  double Inlet_Angle = 2.0 * atan(1.0);


  /// Function that prescribes the hydrostatic pressure field at the outlet
  void hydrostatic_pressure_outlet(const double &time, const Vector<double> &x, const Vector<double> &n,
                                   Vector<double> &traction) {
    traction[0] = ReInvFr * G[1] * (1.0 - x[1]);
    traction[1] = 0.0;
  }

  /// Function that prescribes hydrostatic pressure field at the inlet
  void hydrostatic_pressure_inlet(const double &time, const Vector<double> &x, const Vector<double> &n,
                                  Vector<double> &traction) {
    traction[0] = -ReInvFr * G[1] * (1.0 - x[1]);
    traction[1] = 0.0;
  }

  //end of traction functions

  /// Constitutive law used to determine the mesh deformation
  ConstitutiveLaw *Constitutive_law_pt;

  /// Pseudo-solid Poisson ratio
  double Nu = 0.1;
}

#endif //GLOBALPHYSICALVARIABLES_H
