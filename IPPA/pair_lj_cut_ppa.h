/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.


    Pair-style that continuously performs PPA as a control parameter is varied

    Arguments 
    
        window  integer          Chemical distance to switch off pair interactions
                                 Default 10
        
        alpha   double>0         Exponent for stretching switching/compression process
                                 Default 1
        
        u0      double>=0        Maximal energy at overlap for calculating the inner cutoff
     or beta0   0<=double<=1     Inner cutoff
                                 Default beta=0.817505 corresponding to U0=100
     
        lambda  0<=double<=1     Initial value for mixing parameter.  (before fix adapt)
                                 Default 1.0  (full KG force field)
        
        circular                 Assume chains have circular connectivity
                                 Default false   assumes linear chains.

------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(wca/ppa,PairWCAPPA)

#else

#ifndef LMP_PAIR_WCAPPA
#define LMP_PAIR_WCAPPA

#include "pair.h"

namespace LAMMPS_NS {

class PairWCAPPA : public Pair {
 public:
  PairWCAPPA(class LAMMPS *);
  virtual ~PairWCAPPA();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
//  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);
  void EstimateMoleculeMinMax();

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double *cut_respa;
  
  double nwindow;
  double lambda;
  double alpha,beta,gamma;
  bool circular;
  bool debug;

  int *moleculebeadmax; // maximal tag no in molecule
  int *moleculebeadmin; // minimal tag no in molecule

  virtual void allocate();
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
