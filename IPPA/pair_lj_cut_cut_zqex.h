/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- 

lj/cut with the modification that an inner cut-off is also
defined such that 

    Up(r)=  (r-rc) U'(rc)+ U(rc)  for r<rc
                   U(r)           for r>=rc
                   
   where U=Uljcut(r)
   
Hence the force is

    Fp(r) =   Flj(rc) for r<rc
              Flj(r)  for r>=rc
                     

Reading force field parameters from restart files is not supported!

Syntax:
pair_style lj/cutcut  <globalcutoff>
pair_coeff <i> <j>   <epsilon> <sigma> <rcut_inner> [<rcout_outer>]

where the global cutoff is used if the outer cutoff is not specified.
*/

#ifdef PAIR_CLASS

PairStyle(lj/cutcut,PairLJCutCut)

#else

#ifndef LMP_PAIR_LJ_CUTCUT_H
#define LMP_PAIR_LJ_CUTCUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCut : public Pair {
 public:
  PairLJCutCut(class LAMMPS *);
  virtual ~PairLJCutCut();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut,**cutcut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**lj5,**offset,**offset2;
  double **cutcutsq;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

*/
