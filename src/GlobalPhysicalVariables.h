//
// Created by Oscar Holroyd on 09/01/2025.
//

#ifndef GLOBALPHYSICALVARIABLES_H
#define GLOBALPHYSICALVARIABLES_H

//Finite-Element library routines
#include "generic.h"

using namespace std;

using namespace oomph;

//The global physical variables
namespace Global_Physical_Variables {
  /// The length of the domain in the x direction (streamwise)
  double Lx = 32.0;

  /// The length of the domain in the y direction (spanwise)
  double Ly = 32.0;

  /// Reynolds number, based on the surface velocity of a flat film (Nusselt velocity)
  double Re = 15.0;

  /// The Capillary number
  double Ca = 0.05;

  /// Angle of incline of the slope
  double Theta = M_PI / 4.0;

  /// The Vector direction of gravity (x/y/z, ie streamwise, spanwise, normal)
  Vector<double> G = {2.0, 0.0, -2.0 / tan(Theta)}; // x is streamwise, y is spanwise, z is normal
  // Vector<double> G = {0.0, 2.0, -2.0 / tan(Theta)}; // x is spanwise, y is streamwise, z is normal
  // Vector<double> G = {sqrt(2.0), sqrt(2.0), -2.0 / tan(Theta)}; // x+y is streamwise, x-y is spanwise, z is normal

  /// The product of Reynolds number and inverse Froude number (always 1 in this scaling)
  double ReInvFr = 1.0;
}

#endif //GLOBALPHYSICALVARIABLES_H
