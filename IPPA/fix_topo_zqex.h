/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   Modified from MISC/fix_srp.h

   The fix checks for topology violations for bonds crossing.
   
   Fix exports 4 integers
   
   1:  Instantaneous  (intra and inter) topology violations
   2:  Accumulated    (intra and inter) topology violations
   3:  Instantaneous  (inter only) topology violations
   4:  Accumulated    (inter only) topology violations
   
   Thus most of the time you want option two.   

------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(topo,FixTopo);
// clang-format on
#else

#ifndef LMP_FIX_TOPO_H
#define LMP_FIX_TOPO_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTopo : public Fix {
 public:
  FixTopo(class LAMMPS *, int, char **);
  ~FixTopo() override;
  int setmask() override;
  void init() override;
  void setup_pre_force(int) override;

  void pre_exchange() override;
  void post_run() override;
  void post_force(int) override;

  void min_pre_exchange() override;
  void min_post_force( int ) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_border(int, int *, double *) override;
  int unpack_border(int, int, double *) override;

  int pack_reverse_comm(int n, int first, double *buf) override;
  void unpack_reverse_comm(int n, int *list, double *buf) override;

  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int maxsize_restart() override;
  int size_restart(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  int modify_param(int, char **) override;

// Added by zqex:  
  double compute_vector(int) override;
  double **array;
  bool *flip;
  bool doFlip;
  bool debug;
  
  // All topology violations
  int topologyviolations;
  int topologyviolationsaccumulator;

  // intermolecular:
  int topologyviolationsinter;
  int topologyviolationsinteraccumulator;

 protected:
  int btype;
  int bptype;
  std::string pair_name;
};

}    // namespace LAMMPS_NS

#endif
#endif
