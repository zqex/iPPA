/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Modified from ../MISC/pair_srp.h

    Pair-style that continuously performs PPA as a control parameter is varied

    Arguments 
    
        cutoff  double>0         Cutoff for looking for topology violating bonds.
    
        window  integer          Chemical distance to switch off pair interactions
                                 Default 10
        
        alpha   double>0         Exponent for iPPA switching potential
                                 Default 1
        
        u0      double>=0        Maximal energy at overlap for calculating the inner cutoff
     or beta0   0<=double<=1     Inner cutoff
                                 Default beta=0.817505 corresponding to U0=100
     
        lambda  0<=double<=1     Initial value for mixing parameter.  (before fix adapt modifies it)
                                 Default 1.0  (full KG force field)

        beadmap  <file>          Assigns an integer to each bead. Beads with same integer 
                                 can be contracted. molecule number is used if beadmap is
                                 not specificied but this can lead to loss of entanglement
                                 if used to perform PPA on polymer networks.
        
        circular                 Assume chains have circular connectivity
                                 Default false   assumes linear chains.
   
        flip                     Flips bond pairs that has crossed back to their original position, and flips velocity vectors.
                                 (Experimental feature!)

        debug                    Save topology violations to a per-core debug file.
                                 Default false.
   
------------------------------------------------------------------------- */


#ifdef PAIR_CLASS
// clang-format off
PairStyle(topo,PairTopo);
// clang-format on
#else



#ifndef LMP_PAIR_TOPO_H
#define LMP_PAIR_TOPO_H

#include "pair.h"
#include "atom.h"

namespace LAMMPS_NS {

class PairTopo : public Pair {
 public:
  PairTopo(class LAMMPS *);
  ~PairTopo() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void *extract(const char *, int &) override;
  void EstimateMoleculeMinMax();


 protected:
  inline void onetwoexclude(int *&, int &, int *&, int *&, int **&);
  inline void remapBonds(int &);
  void allocate();
  void getMinDist(double **&, double &, double &, double &, double &, double &, int &, int &, int &,
                  int &);
                  
  double **cut;
  double cut_global;
  int bptype;
  int btype;
  class FixTopo *f_topo;
  int exclude, maxcount;
  int **segment;


// chemical distance window:
  double nwindow;
  double alpha,lambda;
  bool circular;
  bool doFlip;
  bool debug;

  int *moleculebeadmax; // maximal tag no in molecule
  int *moleculebeadmin; // minimal tag no in molecule
    
  // Use map file to specify molecule, index, topology for every bead
  // PPA is applied within molecule, depending on chemical distance
  // based on index differences, and correcting for circular chains
  // if topology is circular.

  bool hasMap;                
  void PPAreadmap(char* str);
  int *molmap,*indexmap;      
  bool *topomap;

  // resolve lookup
  inline int getmol(int i)
          {
              if (hasMap) return molmap[atom->tag[i]];
                     else return atom->molecule[i];
          }

  // resolve lookup
  inline int getindex(int i)
          {
              if (hasMap) return indexmap[atom->tag[i]];
                     else return atom->tag[i];
          }

  // resolve lookup
  inline bool gettopo(int i)
          {
              if (hasMap) return topomap[atom->tag[i]];
                     else return circular;
          }


};

}    // namespace LAMMPS_NS

#endif
#endif
