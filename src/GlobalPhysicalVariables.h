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
  double Re = 15.0;

  /// The product of Reynolds number and inverse Froude number
  double ReInvFr = 1.0;

  /// Angle of incline of the slope (45 degrees)
  double Alpha = M_PI / 4.0;

  /// The Vector direction of gravity, set in main()
  Vector<double> G(2, 0.0);

  /// The Capillary number
  double Ca = 0.05;

  /// Set the wavenumber
  double K = 1;

  /// The length of the domain to fit the desired number of waves
  double Length = 20;
}

#endif //GLOBALPHYSICALVARIABLES_H
