/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
   
Modified by ZQEX:

To perform PPA on KG melts we disable pair interactions between beads of a
certain chemical distance on the same chain (n ~ 10, depending on stiffness),
while retaining the full interaction between beads on different chains.

To invert the PPA in a continuous way, we need a way to continously reintroduce
the interactions and heat up the system such that we can switch to the KG
force field again.

Assuming pair special 1 1 1 and FENE with only spring no WCA.

For x in [0:1] 

x=0 corresponds to PPA with no pair interaction for beads |i-j|<=Nwindow
x=1 corresponds to WCA between all beads.

Progressively WCA is switched on using a force capped potential for bonds from
Nwindow down to the first bond.

Normal pair forces for |i-j|>n(x) with   n(x) = int( (1-x)^alpha*nwindow )

For |i-j|=n(x) we use force capped WCA interaction with an inner cap Rc= 2^(1/6)*sigma*y 
with y=(1-x)*n-n(x). 

For |i-j|<n(x) there are no interactions.

------------------------------------------------------------------------- */

#include "pair_lj_cut_ppa.h"
#include "citeme.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

#include <iostream>
#include <cmath>
#include <cstring>

static const char cite_topo[] =
  "pair topo command: doi:XX\n\n"
  "@Article{Svaneborg2023\n"
  " author = {C. Svaneborg},\n"
  " title = {},\n"
  " journal = {Comp. Phys. Comm.},\n"
  " year =    2023,\n"
  " volume =  xxx,\n"
  " pages =   {xxx}\n"
  "}\n\n";


#include <string.h>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairWCAPPA::PairWCAPPA(LAMMPS *lmp) : Pair(lmp), moleculebeadmax(nullptr), moleculebeadmin(nullptr), debug(false)
{
  writedata = 1;
  single_enable=0;

  cut_global = pow(2.0, 1.0/6.0);   // WCA cutoff. This is fixed. wca/ppa only works with WCA particles.

  nwindow=10;       // Nwindow
  gamma=0.817505;   // Ucap=100
  alpha=1.0;        // No warping
  beta=1.0;         // No warping
  lambda=1.0;       // Starting in KG force field.
  circular = false; // Linear molecules

  if (lmp->citeme) lmp->citeme->add(cite_topo);
}

/* ---------------------------------------------------------------------- */

PairWCAPPA::~PairWCAPPA()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    
    if (moleculebeadmax) memory->destroy(moleculebeadmax);
    if (moleculebeadmin) memory->destroy(moleculebeadmin);
  }
}

/* ---------------------------------------------------------------------- */

void PairWCAPPA::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int *tag = atom->tag;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  tagint *molecule = atom->molecule;

  if (circular && (!moleculebeadmin or !moleculebeadmax)) EstimateMoleculeMinMax();
  if (!list) error->all(FLERR,"pair lj/cut/ppa called without a pair list");

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  double W=pow(1-lambda,alpha)*nwindow+1;
  
  int w=floor( W );
  if (lambda==0) w=nwindow;
  if (lambda==1) w=0;
  
  double xx=pow(W-w,beta);
  double gamma0=gamma+(1-gamma)*xx;                     
   
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq>=cutsq[itype][jtype]) continue;

      double rcuti=gamma0*sigma[itype][jtype]*cut_global; 
      double rcuti2;
      int cdist;
      bool fullLJ;
      bool capLJ;

      if (lambda==1.0) { fullLJ=true; capLJ=false; }             // in this special case use full FF
        else
      if (molecule[i]!=molecule[j]) 
        {
           fullLJ=true; capLJ=false;
        }
      else
        {
              cdist=abs( tag[i]-tag[j] );                 // Chemical distance assuming consequtive bead numbers.
              
              if (circular)  // Also look for chemical distance via end-to-end closure
                 {
                      int cdist2=MIN(tag[i],tag[j])-moleculebeadmin[molecule[i]] + moleculebeadmax[molecule[i]]-MAX(tag[i],tag[j])+1;
                      if (cdist2 < cdist) cdist=cdist2;
                 } 
              
              rcuti2=rcuti*rcuti;
 
              fullLJ =   cdist>w  || ( cdist==w && rsq>=rcuti2 );
              capLJ  =   cdist==w && rsq<rcuti2;
        }

      if (!fullLJ && !capLJ) continue;                  // no interaction.
      if (fullLJ && capLJ)      error->all(FLERR,"pair lj/cut/ppa impossible state");

      double r,lj5,offset2;

      if (fullLJ) { // Calculate WCA interaction
               r2inv = 1.0/rsq;
               r6inv = r2inv*r2inv*r2inv;
               forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
               fpair = factor_lj*forcelj*r2inv;

               if (eflag) {
                   evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) - offset[itype][jtype];
                   evdwl *= factor_lj;
                  }
            }
      if (capLJ)
            {         // Calculate forcecap.
               r=sqrt(rsq);
               double ratio=sigma[itype][jtype] / rcuti;
               double ratio2=ratio*ratio;
               double ratio4=ratio2*ratio2;
               double ratio6=ratio4*ratio2;
               double ratio12=ratio6*ratio6;

               lj5= 24*epsilon[itype][jtype]*( -ratio6 + 2*ratio12 )/rcuti; 

               fpair = factor_lj*lj5/r;

                if (eflag) {
                   double offset2=4.0 * epsilon[itype][jtype] * (ratio12- ratio6)-offset[itype][jtype];
                   evdwl = factor_lj*(  (rcuti-r)*lj5+offset2  );
               }
           }
       
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}



void PairWCAPPA::EstimateMoleculeMinMax()
/*
  Find maximal & minimal bead numbers in all molecules. This is required for calculating chemical distances
  for circular molecules.
*/
{
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
   allocate all arrays
------------------------------------------------------------------------- */

void PairWCAPPA::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairWCAPPA::settings(int narg, char **arg)
{
/*
   utils::logmesg(lmp,"Pair wca/ppa:  narg={}\n",narg);  
   for (int i=0;i<narg;i++)
       utils::logmesg(lmp,"Pair wca/ppa:  arg[{}]={}\n",i,arg[i]);  
*/

  int iarg=0;
  while (iarg < narg) {
      if (strcmp(arg[iarg],"window") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style wca/ppa command. Missing window argument");
      nwindow    = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (nwindow<=0) error->all(FLERR,"Illegal pair_style wca/ppa command. nwindow>0!");
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"alpha") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style wca/ppa command. Missing alpha argument");
      alpha      = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (alpha<=0) error->all(FLERR,"Illegal pair_style wca/ppa command. alpha>0!");
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"beta") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style wca/ppa command. Missing beta argument");
      beta      = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (beta<=0) error->all(FLERR,"Illegal pair_style wca/ppa command. beta>0!");
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"u0") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style wca/ppa command. Missing U0 argument");
      double u0      = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (u0<0) error->all(FLERR,"Illegal pair_style wca/ppa command. u0=>0!");
      gamma=pow(13.0/(7.0+sqrt(36.0+13.0*u0)),1.0/6.0);
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"gamma") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style wca/ppa command. Missing gamma argument");
      gamma  = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (gamma<0) error->all(FLERR,"Illegal pair_style wca/ppa command. gamma<0!");
      if (gamma>1) error->all(FLERR,"Illegal pair_style wca/ppa command. gamma>1!");
      iarg += 2; }
    else 
      if (strcmp(arg[iarg],"lambda") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal pair_style wca/ppa command. Missing lambda argument");
      lambda  = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (lambda<0) error->all(FLERR,"Illegal pair_style wca/ppa command. lambda<0!");
      if (lambda>1) error->all(FLERR,"Illegal pair_style wca/ppa command. lambda>1!");
      iarg += 2; }
     else 
      if (strcmp(arg[iarg],"circular") == 0) {
      circular=true;      
      iarg += 1; }        
     else 
      if (strcmp(arg[iarg],"debug") == 0) {
      debug=true;      
      iarg += 1; }        
     else error->all(FLERR,fmt::format("Illegal pair wca/ppa command. Unknown argument {}",arg[iarg]) );
  }
  
  double U0=1.0+13.0/pow(gamma,12)-14.0/pow(gamma,6);
  if (comm->me==0)
    {
      utils::logmesg(lmp,fmt::format("Pair WCA/ppa: nwindow={}\n",nwindow)  );  
      utils::logmesg(lmp,fmt::format("Pair WCA/ppa: alpha={}\n",alpha)  );  
      utils::logmesg(lmp,fmt::format("Pair WCA/ppa: beta={}\n",beta)  );  
      utils::logmesg(lmp,fmt::format("Pair WCA/ppa: U0={}\n",U0)  );  
      utils::logmesg(lmp,fmt::format("Pair WCA/ppa: lambda={}\n",lambda)  );
      if (circular)   utils::logmesg(lmp, "Pair WCA/ppa: circular chains\n" );
                else  utils::logmesg(lmp, "Pair WCA/ppa: linear chains\n" );
    }
    
  // reset cutoffs that have been explicitly set
  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairWCAPPA::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR,arg[4],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
    
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairWCAPPA::init_style()
{
    if (atom->molecule == NULL)
    error->all(FLERR,
               "Must use atom style with molecule IDs with pair lj/cut/ppa");

  if (force->special_lj[1] != 1.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0)
    error->all(FLERR,"pair lj/cut/ppa requires special_bonds = 1,1,1");

  // request regular neighbour list
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairWCAPPA::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  sigma[j][i] = sigma[i][j];
  epsilon[j][i] = epsilon[i][j];

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWCAPPA::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWCAPPA::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairWCAPPA::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairWCAPPA::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairWCAPPA::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairWCAPPA::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairWCAPPA::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  dim=0;
  if (strcmp(str,"lambda") == 0) return (void *)&lambda;
  return nullptr;
}
