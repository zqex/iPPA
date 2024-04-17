/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------

  poly06 is a 6th order taylor expansion of a bond around the equilibrium distance r0.

  U(r)= k0 + k1(r-r0) + k2(r-r0)^2 + k3(r-r0)^3 + k4(r-r0)^4  + k5(r-r0)^5  + k6(r-r0)^6
  
  bond_style poly06
  bond_coeff *  <r0>  <k0> <k1> <k2> <k3> <k4> <k5> <k6>

 */

#ifdef BOND_CLASS
// clang-format off
BondStyle(poly06,BondPoly06);
// clang-format on
#else

#ifndef LMP_BOND_POLY06_H
#define LMP_BOND_POLY06_H

#include "bond.h"

namespace LAMMPS_NS {

class BondPoly06 : public Bond {
 public:
  BondPoly06(class LAMMPS *);
  ~BondPoly06() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *r0,*k0,*k1,*k2,*k3,*k4,*k5,*k6;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
