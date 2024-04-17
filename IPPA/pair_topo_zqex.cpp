// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Timothy Sirk (ARL), Pieter in't Veld (BASF)

This pair style srp command calculates a segmental repulsive force
between bonds. This is useful for preventing the crossing of bonds if
soft non-bonded potentials are used, such as DPD polymer chains.

See the doc page for pair_style srp command for usage instructions.

There is an example script for this package in examples/PACKAGES/srp.

Please contact Timothy Sirk for questions (tim.sirk@us.army.mil).
-------------------------------------------------------------------------

  Modified from ../MISC/pair_srp.cpp
  
  The SRP code bond has been removed, but this code is combined with fix topo to
  count topology violations.
  
*/

#include "pair_topo_zqex.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_topo_zqex.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#include "thermo.h"

#include "update.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

using namespace LAMMPS_NS;

#define SMALL 1.0e-10
#define BIG 1e10
#define ONETWOBIT 0x40000000

static const char cite_srp[] =
  "pair srp command: doi:10.1063/1.3698476\n\n"
  "@Article{Sirk2012\n"
  " author = {T. W. Sirk and Y. R. Sliozberg and J. K. Brennan and M. Lisal and J. W. Andzelm},\n"
  " title = {An Enhanced Entangled Polymer Model for Dissipative Particle Dynamics},\n"
  " journal = {J.~Chem.\\ Phys.},\n"
  " year =    2012,\n"
  " volume =  136,\n"
  " pages =   {134903}\n"
  "}\n\n";

static const char cite_topo[] =
  "pair topo command: doi:XX\n\n"
  "@Article{Svaneborg2024\n"
  " author = {C. Svaneborg},\n"
  " title = {},\n"
  " journal = {Comput. Phys. Commun.},\n"
  " year =    2024,\n"
  " volume =  xxx,\n"
  " pages =   {xxx}\n"
  "}\n\n";

/* ----------------------------------------------------------------------
 set size of pair comms in constructor
 ---------------------------------------------------------------------- */

PairTopo::PairTopo(LAMMPS *lmp) : Pair(lmp), moleculebeadmax(nullptr), moleculebeadmin(nullptr)
{
  writedata = 1;
  single_enable = 0;

  if (lmp->citeme)
    {
       lmp->citeme->add(cite_srp);
       lmp->citeme->add(cite_topo);
    }

  nextra = 1;
  segment = nullptr;
  circular=false;
  doFlip=false;
  debug=false;

  // create fix SRP instance here with unique fix id
  // similar to granular pair styles with history,
  //   this should be early enough that FixSRP::pre_exchange()
  //   will be invoked before other fixes that migrate atoms
  //   this is checked for in FixSRP

  if (modify->nfix && modify->get_fix_by_id("topo")!=nullptr)
    {
       f_topo = dynamic_cast<FixTopo*>(modify->get_fix_by_id("topo"));
       if (comm->me==0) 
            utils::logmesg(lmp,"PairTopo: Reusing old fix topo all topo\n");
    }
   else
    {   
        f_topo = dynamic_cast<FixTopo*>(modify->add_fix(fmt::format("topo all topo")));
        if (comm->me==0) 
            utils::logmesg(lmp,"PairTopo: Creating new: fix topo all topo\n");
    }

}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairTopo::allocate()
{
    allocated = 1;
    // particles of bptype inserted by fix srp
    // bptype is the highest numbered atom type
    int n = bptype;
    memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
    memory->create(cut, n + 1, n + 1, "pair:cut");

    // setflag for atom types
    memory->create(setflag,n+1,n+1,"pair:setflag");
    for (int i = 1; i <= n; i++)
        for (int j = i; j <= n; j++)
            setflag[i][j] = 0;

    maxcount = 0;
}

/* ----------------------------------------------------------------------
 free
 ------------------------------------------------------------------------- */

PairTopo::~PairTopo()
{

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(segment);
  }
  if (moleculebeadmax) memory->destroy(moleculebeadmax); moleculebeadmax=nullptr;
  if (moleculebeadmin) memory->destroy(moleculebeadmin); moleculebeadmin=nullptr;

  // check nfix in case all fixes have already been deleted
  if (modify->nfix && modify->get_fix_by_id(f_topo->id)!=nullptr)
     { 
       modify->delete_fix(f_topo->id);
       if (comm->me==0) 
            utils::logmesg(lmp,"PairTopo: Destroying fix topo\n");
     }
}

/* ----------------------------------------------------------------------
 compute bond-bond repulsions
 ------------------------------------------------------------------------- */

void PairTopo::compute(int eflag, int vflag)
{
    // setup energy and virial
    ev_init(eflag, vflag);

    double **x = atom->x;
    double **v = atom->v;
    int nlocal = atom->nlocal;
    int nall = nlocal + atom->nghost;
    int i0, i1, j0, j1;
    int i,j,ii,jj,inum,jnum;
    double dijsq, dij;

    int *ilist,*jlist,*numneigh,**firstneigh;
    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    double dx,dy,dz,ti,tj;
    double **xlast= f_topo->array;
    int *tag = atom->tag;
    tagint *molecule = atom->molecule;
    
    // mapping global to local for atoms inside bond particles
    // exclude 1-2 neighs if requested
    if (neighbor->ago == 0) {
      remapBonds(nall);
      if (exclude) onetwoexclude(ilist, inum, jlist, numneigh, firstneigh);
    }

  // Current chemical window w
  double W=pow(1-lambda,alpha)*nwindow+1;  
  int w=floor( W );
  if (lambda==0) w=nwindow;
  if (lambda==1) w=0;

  if (circular && (!moleculebeadmin || !moleculebeadmax))
     {
        if (debug) 
           utils::logmesg(lmp,fmt::format("PairTopo: compute  circular={}  moleculebeadmin={} moleculebeadmax={}\n",circular,(void*)moleculebeadmin,(void*)moleculebeadmax));

        EstimateMoleculeMinMax();
     }

  // counter for topology violations at this time step.
  f_topo->topologyviolations=0;
  f_topo->topologyviolationsinter=0;

  // this pair style only used with hybrid
  // due to exclusions
  // each atom i is type bptype
  // each neigh j is type bptype

  int mol0,mol1;

  for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];
      jnum = numneigh[i];
      i0 = segment[i][0];
      j0 = segment[i][1];
      if (atom->type[i] != bptype) error->all(FLERR,"Illegal bond particle type in topology loop");

//      utils::logmesg(lmp,fmt::format("topo {}:{}: i={} tag={} at {} {} {}   {}\n",update->ntimestep,comm->me,  i,atom->tag[i],x[i][0],x[i][1],x[i][2], atom->type[i]));

      for (jj = 0; jj < jnum; jj++) {

        jlist = firstneigh[i];
        j = jlist[jj];
        j &= NEIGHMASK;

        if (atom->type[j] != bptype) error->all(FLERR,"Illegal bond particle type in topology loop");

        i1 = segment[j][0];
        j1 = segment[j][1];

        // Find molecule of bond. This is only defined if both beads in the bond belong to the same molecule.
        mol0=-1;
        if (molecule[i0]==molecule[j0]) mol0=molecule[j0];

        mol1=-1;
        if (molecule[i1]==molecule[j1]) mol1=molecule[j1];

        // Topology check, only worry about chemical distance, if bonds can be assigned to molecules, and
        //  they belong to the same molecule. Then calculate chemical distance and compare to current chemical window
                    
        if (mol0>-1 && mol1>-1 && mol0 == mol1)         // Intramolecular bond pair.
           {
              // Chemical distance, since we do not know the order of beads in the two bond segments try all options.
              int cdist=abs(tag[i0]-tag[i1]);           
              if (abs(tag[i0]-tag[j1])<cdist) cdist=abs(tag[i0]-tag[j1]);
              if (abs(tag[j0]-tag[i1])<cdist) cdist=abs(tag[j0]-tag[i1]);
              if (abs(tag[j0]-tag[j1])<cdist) cdist=abs(tag[j0]-tag[j1]);

              // If circular we also have to look at distances across the starting and ending beads.
              if (circular && cdist>0)                  // then look for chemical distance via end-to-end closure,
                 {                                      // if cdist==0 when one of the i0,j0 i1,j1 beads are the same, so two neighboring bonds.
                      int mintag=MIN( MIN(tag[i0],tag[j0]), MIN(tag[i1],tag[j1]) );
                      int maxtag=MAX( MAX(tag[i0],tag[j0]), MAX(tag[i1],tag[j1]) );
                                       
                      //           length of head                 length of tail                and closure bond
                      int cdist2=  mintag-moleculebeadmin[mol0]  +moleculebeadmax[mol0]-maxtag  +1;
                      if (cdist2 < cdist) cdist=cdist2;
                 } 

              if (cdist<=w) continue;            // if chemical distance smaller than current window, do not check topology.
           }

  
        // If we flip bond that cross, do not compare topology violation if bond was flipped in last step.     
        if (doFlip && x[i0][0]==xlast[i0][0] && x[j0][0]==xlast[j0][0] && x[i1][0]==xlast[i1][0] && x[j1][0]==xlast[j1][0]   &&
                      x[i0][1]==xlast[i0][1] && x[j0][1]==xlast[j0][1] && x[i1][1]==xlast[i1][1] && x[j1][1]==xlast[j1][1]   &&
                      x[i0][2]==xlast[i0][2] && x[j0][2]==xlast[j0][2] && x[i1][2]==xlast[i1][2] && x[j1][2]==xlast[j1][2] ) continue;


        // This gets the vector dx,dy,dz in the current config connecting the two segments defined by
        // bead pair i0,j0  and bead pair i1,j1
        getMinDist(x, dx, dy, dz, ti, tj, i0, j0, i1, j1);
        if (ti<-0.5 or ti>0.5) continue; // Skip bond crossing check, if min distances occurs outside segment.
        if (tj<-0.5 or tj>0.5) continue; // Skip bond crossing check, if min distances occurs outside segment.

        dijsq = dx*dx + dy*dy + dz*dz;
        if (dijsq >= cutsq[bptype][bptype]) continue;
      
      // Identify closest images of xlast in relation to present coordinates.
        domain->remap_near(xlast[i0],x[i0]);
        domain->remap_near(xlast[i1],x[i1]);
        domain->remap_near(xlast[j0],x[j0]);
        domain->remap_near(xlast[j1],x[j1]);

        double dxlast,dylast,dzlast;
        getMinDist(xlast, dxlast, dylast, dzlast, ti, tj, i0, j0, i1, j1);
        if (ti<-0.5 or ti>0.5) continue;
        if (tj<-0.5 or tj>0.5) continue;

        // We have a bond violation event.
        if (dx*dxlast+dy*dylast+dz*dzlast < 0)
            {
               if (debug)
                  {
                     std::string fnam = fmt::format("DEBUG_topo.me{}.txt",comm->me);
                     std::ofstream fo(fnam.c_str(),std::ios_base::app);
                     fo << fmt::format("{} VIOLATION count={}  lambda={}   {}-{} and {}-{}  on {},{}   at {},{},{}    dot={}\n", update->ntimestep, f_topo->topologyviolationsinteraccumulator, lambda, tag[i0],tag[j0],  tag[i1],tag[j1], 
                     mol0,mol1 , 0.25*(x[i0][0]+x[j0][0]+x[i1][0]+x[j1][0])  
                               , 0.25*(x[i0][1]+x[j0][1]+x[i1][1]+x[j1][1])
                               , 0.25*(x[i0][2]+x[j0][2]+x[i1][2]+x[j1][2])
                               ,dx*dxlast+dy*dylast+dz*dzlast
                                );
                                
                     for (int i=tag[i0]-5 ; i<=tag[j0]+5; i++)
                       {
                           int iii = atom->map(i);
                           if (iii==-1) continue;
                           if (molecule[iii]!=molecule[i0]) continue;

                           fo << "\t\t" << i
                              << "     " << xlast[iii][0] << " " << xlast[iii][1] << " " << xlast[iii][2]
                              << "     " << x    [iii][0] << " " << x    [iii][1] << " " << x    [iii][2]
                              << "\n";
                       }

                     for (int i=tag[i1]-5 ; i<=tag[j1]+5; i++)
                       {
                           int iii = atom->map(i);
                           if (iii==-1) continue;
                           if (molecule[iii]!=molecule[i1]) continue;

                           fo << "\t\t" << i
                              << "     " << xlast[iii][0] << " " << xlast[iii][1] << " " << xlast[iii][2]
                              << "     " << x    [iii][0] << " " << x    [iii][1] << " " << x    [iii][2]
                              << "\n";
                       }

                                
                     fo.close();
                  }

               // Sample all topology violations
               f_topo->topologyviolations++;
               f_topo->topologyviolationsaccumulator++;                

               if (!(mol0>-1 && mol1>-1 && mol0 == mol1) )  // intermolecular bond breakage
                  {
                     f_topo->topologyviolationsinter++;
                     f_topo->topologyviolationsinteraccumulator++;
                  }
                                                                      
               // Should be performed in a parallel manner!!
               // Flips involving ghost particles, should be communicated back to their parent particles.
               if (doFlip)
                 {
                    f_topo->flip[i0]=true;
                    f_topo->flip[i1]=true;
                    f_topo->flip[j0]=true;
                    f_topo->flip[j1]=true;
                 }
            }
      }
  }
}


void PairTopo::EstimateMoleculeMinMax()
/*
  Find maximal & minimal bead numbers in all molecules. This is required for calculating chemical distances
  for circular molecules.
*/
{
  if (debug) 
            utils::logmesg(lmp,fmt::format("PairTopo: Calling estimate molecule min/max me={}",comm->me));

  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
   
  // Fix maximal molecule number globally:
  int moleculemax=0;
  int localmoleculemax=0;
  for (int i = 0; i < nlocal; i++)
       {
         if (molecule[i]>localmoleculemax) localmoleculemax=molecule[i];
       }
  MPI_Allreduce(&localmoleculemax,&moleculemax,1,MPI_INT,MPI_MAX,world);
  moleculemax++;  // One larger than the largest molecule, makes space if there is a zero molecule.
   
  // For each molecule find the largest and smallest ID number and store them.
  // We assume this are the ends of a linear molecule

  int *localmoleculebeadmax, *localmoleculebeadmin;
  memory->create(localmoleculebeadmax,moleculemax,"pair_lj_cut_ppa init localmoleculebeadmax");
  memory->create(localmoleculebeadmin,moleculemax,"pair_lj_cut_ppa init localmoleculebeadmin");

  if (moleculebeadmax) memory->destroy(moleculebeadmax);
  if (moleculebeadmin) memory->destroy(moleculebeadmin);
  memory->create(moleculebeadmax,moleculemax,"pair_lj_cut_ppa init moleculebeadmax");
  memory->create(moleculebeadmin,moleculemax,"pair_lj_cut_ppa init moleculebeadmin");
      
  for (int i=0;i<moleculemax;i++) localmoleculebeadmax[i]=0;
  for (int i=0;i<moleculemax;i++) localmoleculebeadmin[i]=atom->map_tag_max;

  for (int i = 0; i < nlocal; i++)
     {
        int m=molecule[i];
        int t=tag[i];
        if (t<localmoleculebeadmin[m]) localmoleculebeadmin[m]=t;
        if (t>localmoleculebeadmax[m]) localmoleculebeadmax[m]=t;
     }
        
  MPI_Allreduce(localmoleculebeadmin,moleculebeadmin,moleculemax,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(localmoleculebeadmax,moleculebeadmax,moleculemax,MPI_INT,MPI_MAX,world);
    
  memory->destroy(localmoleculebeadmin);
  memory->destroy(localmoleculebeadmax);

  if (comm->me==0 && debug)
     for (int m = 1; m < moleculemax; m++)
         utils::logmesg(lmp,fmt::format("Pair wca/ppa: Molecule  {}  from {} to {}\n",m,moleculebeadmin[m], moleculebeadmax[m] ));  
} 


/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairTopo::settings(int narg, char **arg)
{
/*   
   utils::logmesg(lmp,"Pair topo:  narg={}\n",narg);  
   for (int i=0;i<narg;i++)
       utils::logmesg(lmp,"Pair topo:  arg[{}]={}\n",i,arg[i]);  
*/
  int iarg=0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style topo command.  Missing cutoff argument");
      cut_global = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2; }
    else 
    if (strcmp(arg[iarg],"bondtype") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style topo command.  Missing bondtype argument");
      
      if (strcmp(arg[iarg+1],"*") == 0) { btype = 0; }
        else {
              btype = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
              if ((btype > atom->nbondtypes) || (btype <= 0))
                      error->all(FLERR,"Illegal pair_style command");
             }
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"window") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style topo command. Missing window argument");
      nwindow    = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nwindow<=0) error->all(FLERR,"Illegal pair_style topo command. nwindow>0!");
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"alpha") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style topo command. Missing alpha argument");
      alpha      = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (alpha<=0) error->all(FLERR,"Illegal pair_style topo command. alpha>0!");
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"lambda") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style topo command. Missing lambda argument");
      lambda  = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (lambda<0) error->all(FLERR,"Illegal pair_style topo command. lambda<0!");
      if (lambda>1) error->all(FLERR,"Illegal pair_style topo command. lambda>1!");
      iarg += 2; }
     else 
      if (strcmp(arg[iarg],"circular") == 0) {
      circular=true;      
      iarg += 1; }        
     else 
      if (strcmp(arg[iarg],"flip") == 0) {
      doFlip=true;      
      iarg += 1; }        
     else 
      if (strcmp(arg[iarg],"debug") == 0) {
      debug=true;      
      iarg += 1; }        
     else error->all(FLERR,fmt::format("Illegal pair topo command. Unknown argument {}",arg[iarg]) );
  }

  f_topo->doFlip=doFlip;
  f_topo->debug=debug;

  // Current chemical window w
  double W=pow(1-lambda,alpha)*nwindow+1;  
  int w=floor( W );
  if (lambda==0) w=nwindow;
  if (lambda==1) w=0;

  if (comm->me==0)
  {
      utils::logmesg(lmp,fmt::format("Pair topo: nwindow={}\n",nwindow)  );  
      utils::logmesg(lmp,fmt::format("Pair topo: alpha={}\n",alpha)  );  
      utils::logmesg(lmp,fmt::format("Pair topo: lambda={}\n",lambda)  );
      if (circular)   utils::logmesg(lmp, "Pair topo: circular chains\n" );
                else  utils::logmesg(lmp, "Pair topo: linear chains\n" );

      if (doFlip  )   utils::logmesg(lmp, "Pair topo: doFlip enabled\n" );
                else  utils::logmesg(lmp, "Pair topo: doFlip disabled\n" );

      if (debug  )    utils::logmesg(lmp, "Pair topo: debug true\n" );
                else  utils::logmesg(lmp, "Pair topo: debug false\n" );

      utils::logmesg(lmp,fmt::format("Initial window size={}\n",w)  );
  }

  // use last atom type by default for bond particles
  bptype = atom->ntypes;
  exclude = 1;

  // reset cutoffs if explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= bptype; i++)
      for (j = i; j <= bptype; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
 set coeffs
 ------------------------------------------------------------------------- */

void PairTopo::coeff(int narg, char **arg)
{
  if (narg > 3)
    error->all(FLERR,"PairTopo: Incorrect args for pair coeff");
  if (!allocated) allocate();

  // set ij bond-bond cutoffs
  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR,arg[0], 1, bptype, ilo, ihi, error);
  utils::bounds(FLERR,arg[1], 1, bptype, jlo, jhi, error);

  double cut_one = cut_global;
  if (narg == 1)  cut_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      cutsq[i][j] = cut_one * cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->warning(FLERR,"PairTopo: No pair coefficients were set");
}

/* ----------------------------------------------------------------------
 init specific to this pair style
 ------------------------------------------------------------------------- */

void PairTopo::init_style()
{
  if (!force->newton_pair)
    error->all(FLERR,"Pair topo requires newton pair on");

  // verify that fix SRP is still defined and has not been changed.

  if (strcmp(f_topo->style,"topo") != 0)
    error->all(FLERR,"Fix topo has been changed unexpectedly");

  if (comm->me == 0)
    utils::logmesg(lmp,"Using type {} for bond particles\n",bptype);

  // set bond and bond particle types in fix srp
  // bonds of this type will be represented by bond particles
  // if bond type is 0, then all bonds have bond particles
  // btype = bond type

  char c0[20];
  char* arg0[2];
  sprintf(c0, "%d", btype);
  arg0[0] = (char *) "btype";
  arg0[1] = c0;
  f_topo->modify_params(2, arg0);

  // bptype = bond particle type
  sprintf(c0, "%d", bptype);
  arg0[0] = (char *) "bptype";
  arg0[1] = c0;
  f_topo->modify_params(2, arg0);

  // bond particles do not contribute to energy or virial
  // bond particles do not belong to group all
  // but thermo normalization is by nall
  // therefore should turn off normalization
  char *arg1[2];
  arg1[0] = (char *) "norm";
  arg1[1] = (char *) "no";
  output->thermo->modify_params(2, arg1);
  if (comm->me == 0) error->message(FLERR,"Thermo normalization turned off by pair topo");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairTopo::init_one(int i, int j)
{
 if (setflag[i][j] == 0) error->all(FLERR,"PairTopo: All pair coeffs are not set");

  cut[j][i] = cut[i][j];
  return cut[i][j];
}

/* ----------------------------------------------------------------------
 find min distance for bonds i0/j0 and i1/j1
 ------------------------------------------------------------------------- */
inline void PairTopo::getMinDist(double** &x, double &dx, double &dy, double &dz, double &ti, double &tj, int &i0, int &j0, int &i1, int &j1)
{
    // move these outside the loop
    double diffx0, diffy0, diffz0, diffx1, diffy1, diffz1, dPx, dPy, dPz, RiRi, RiRj, RjRj;
    double denom, termx0, termy0, termz0, num0, termx1, termy1, termz1, num1;

    // compute midpt dist from 1st atom, 1st bond
    diffx0 = x[j0][0] - x[i0][0]; // x,y,z from bond 0
    diffy0 = x[j0][1] - x[i0][1];
    diffz0 = x[j0][2] - x[i0][2];

    // compute midpt dist from 1st atom, 2nd bond
    diffx1 = x[j1][0] - x[i1][0];
    diffy1 = x[j1][1] - x[i1][1];
    diffz1 = x[j1][2] - x[i1][2];

    // midpoint distance
    dPx = 0.5*(diffx0-diffx1) + x[i0][0]-x[i1][0];
    dPy = 0.5*(diffy0-diffy1) + x[i0][1]-x[i1][1];
    dPz = 0.5*(diffz0-diffz1) + x[i0][2]-x[i1][2];

    // Ri^2 Rj^2
    RiRi = diffx0*diffx0 + diffy0*diffy0 + diffz0*diffz0;
    RiRj = diffx0*diffx1 + diffy0*diffy1 + diffz0*diffz1;
    RjRj = diffx1*diffx1 + diffy1*diffy1 + diffz1*diffz1;
    denom = RiRj*RiRj - RiRi*RjRj;

    // handle case of parallel lines
    // reduce to midpt distance
    if (fabs(denom) < SMALL) {
        if (denom < 0) denom = -BIG;
        else denom = BIG;
    }

    // calc ti
    termx0 = RiRj*diffx1 - RjRj*diffx0;
    termy0 = RiRj*diffy1 - RjRj*diffy0;
    termz0 = RiRj*diffz1 - RjRj*diffz0;
    num0 = dPx*termx0 + dPy*termy0 + dPz*termz0;
    ti = num0 / denom;

// Changed by zqex:
//    if (ti > 0.5) ti = 0.5;
//    if (ti < -0.5) ti = -0.5;

    // calc tj
    termx1 = RiRj*diffx0 - RiRi*diffx1;
    termy1 = RiRj*diffy0 - RiRi*diffy1;
    termz1 = RiRj*diffz0 - RiRi*diffz1;
    num1 = dPx*termx1 + dPy*termy1 + dPz*termz1;
    tj = -num1/ denom;

// Changed by zqex:
//    if (tj > 0.5)  tj = 0.5;
//    if (tj < -0.5) tj = -0.5;

    // min dist
    dx = dPx - ti*diffx0 + tj*diffx1;
    dy = dPy - ti*diffy0 + tj*diffy1;
    dz = dPz - ti*diffz0 + tj*diffz1;
}

/* --------------------------------------------------------
map global id of atoms in stored by each bond particle
 ------------------------------------------------------- */
inline void PairTopo::remapBonds(int &nall)
{
  if (nall > maxcount) {
    memory->grow(segment, nall, 2, "pair:segment");
    maxcount = nall;
  }

  // loop over all bond particles
  // each bond paricle holds two bond atoms
  // map global ids of bond atoms to local ids
  // might not be able to map both bond atoms of j, if j is outside neighcut
  // these are not on neighlist, so are not used
  int tmp;
  double **srp = f_topo->array_atom;

    for (int i = 0; i < nall; i++) {
      if (atom->type[i] == bptype) {
        // tmp is local id
        // tmp == -1 is ok
        tmp = atom->map((int)srp[i][0]);
        segment[i][0] = domain->closest_image(i,tmp);
        // repeat with other id
        tmp = atom->map((int)srp[i][1]);
        segment[i][1] = domain->closest_image(i,tmp);
      }
    }
}

/* --------------------------------------------------------
add exclusions for 1-2 neighs, if requested
more complex exclusions or scaling probably not needed
 ------------------------------------------------------- */
inline void PairTopo::onetwoexclude(int* &ilist, int &inum, int* &jlist, int* &numneigh, int** &firstneigh)
{
    int i0, i1, j0, j1;
    int i,j,ii,jj,jnum;

    // encode neighs with exclusions
    // only need 1-2 info for normal uses of srp
    // add 1-3, etc later if ever needed

    for (ii = 0; ii < inum; ii++) {

      i = ilist[ii];
      jnum = numneigh[i];
      // two atoms inside bond particle
      i0 = segment[i][0];
      j0 = segment[i][1];

      for (jj = 0; jj < jnum; jj++) {

        jlist = firstneigh[i];
        j = jlist[jj];
        j &= NEIGHMASK;
        //two atoms inside bond particle
        i1 = segment[j][0];
        j1 = segment[j][1];

        // check for a 1-2 neigh
        if (i0 == i1 || i0 == j1 || i1 == j0 || j0 == j1) {
          j |= ONETWOBIT;
          jlist[jj] = j;
        }
      }
    }
}

/* ----------------------------------------------------------------------
proc 0 writes to data file
------------------------------------------------------------------------- */

void PairTopo::write_data(FILE *fp)
{
}

/* ----------------------------------------------------------------------
proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairTopo::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g\n",i,j,cut[i][j]);
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTopo::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTopo::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}
/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairTopo::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&bptype,sizeof(int),1,fp);
  fwrite(&btype,sizeof(int),1,fp);
  fwrite(&exclude,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairTopo::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&bptype,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&btype,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&exclude,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
}

void *PairTopo::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"lambda") == 0) return (void *)&lambda;
  return nullptr;
}
