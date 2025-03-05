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
namespace Global_Variables {
/// The length of the domain in the x direction (streamwise)
double Lx = 32.0;

/// The length of the domain in the y direction (spanwise)
double Ly = 16.0;

/// Reynolds number, based on the surface velocity of a flat film (Nusselt velocity)
double Re = 15.0;

/// The Capillary number
double Ca = 0.05;

/// Angle of incline of the slope
double Theta = M_PI / 4.0;

/// Time to run before the control is turned on
double tburn = 200.0;

/// Time step during burn
double dtburn = 0.1;

/// Time to run with control turned on
double tcontrol = 100.0;

/// Time step during burn
double dtcontrol = 0.1;

/// The number of elements in the x direction
unsigned nx = 40;

/// The number of elements in the y direction
unsigned ny = 20;

/// The number of elements in the z direction
unsigned nz = 2;

/// The number of points in th control system in the x direction
unsigned nx_control = 40;

/// The number of points in th control system in the y direction
unsigned ny_control = 20;

/// The number of actuators in the control system
unsigned m_control = 200;

/// The number of observers in the control system
unsigned p_control = 200;

/// Whether to use the MUMPS solver
bool use_mumps = true;

/// The Vector direction of gravity (x/y/z, ie streamwise, spanwise, normal)
Vector<double> G = {2.0, 0.0, -2.0 / tan(Theta)}; // x is streamwise, y is spanwise, z is normal

/// The product of Reynolds number and inverse Froude number (always 1 in this scaling)
double ReInvFr = 1.0;
}

#endif //GLOBALPHYSICALVARIABLES_H
